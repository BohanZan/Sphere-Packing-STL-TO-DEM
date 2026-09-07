# SpherePacking — Windows C++ / Intel oneAPI

将当前工作目录中的 MATLAB 实现移植为独立 C++20 库和命令行程序。原有 `.m` 文件保持原样。默认构建使用 Intel oneAPI DPC++/C++ Compiler 2025.2，CPU x64、Release 优化和严格浮点运算。

## 编译与运行（PowerShell）

```powershell
Set-Location 'C:\Users\Lenovo\Desktop\Study\AsteroidResearch\DEMBody-move\DEMTools\SpherePacking-cpp'

# 初始化 x64 编译环境、编译并运行回归测试
.\cpp\scripts\Build.ps1

# 当前 testRun.m 参数：请求 22000 个、半径 0.3、重力沿 -X
.\build\release\sphere_packing.exe --preset buddha --report .\results\Great_Budda_cpp\full_report.json
```

构建脚本自动定位 Visual Studio C++ 工具和 oneAPI，优先选择本机 v143/14.44 工具集。运行库会复制到 `.exe` 旁边，运行不需要先进入 oneAPI 命令行。已专门验证普通 PowerShell 以及仅含 Windows 系统路径的环境。

快速运行真实 STL：

```powershell
.\build\release\sphere_packing.exe `
  --model .\inputs\greatBudda\greatBudda.stl `
  --radius 0.625 --count 50 --attempts 60 `
  --sweeps 200 --shake-sweeps 5 --refill-passes 6 `
  --compression-tolerance 1e-7 --gravity 0 0 -1 --frame world `
  --output .\results\cpp_quick --prefix packing `
  --report .\results\cpp_quick\report.json

# 查看所有选项
.\build\release\sphere_packing.exe --help

# 内存访问检查版本；脚本同时部署 AddressSanitizer DLL
.\cpp\scripts\Build.ps1 -Sanitize
.\build\asan\test_packing.exe
```

本机依赖：Intel oneAPI 2025.2、Visual Studio C++ 工具、CMake >= 3.24、Ninja。运行主程序不依赖 MATLAB、Python、网络或第三方 STL 库。MATLAB/NumPy/SciPy 仅用于交叉验证与基准测试。

## 参数与原逻辑

保留有序半径序列、三次空初始尝试结束批次、仅新批次移动、压实与独立侧向扰动、累计空补填遍历上限、任意重力方向、按模型尺度缩放的容差、容量不足时返回有效部分结果。

`--preset buddha` 对应当前 `testRun.m` 的物理和迭代参数；额外固定 C++ 原生种子为 42，输出到 `results/Great_Budda_cpp`。预设后的参数可以覆盖默认值。

不均匀粒径使用无表头的文本/CSV 文件，按顺序给出正半径：

```powershell
.\build\release\sphere_packing.exe --model .\inputs\greatBudda\greatBudda.stl `
  --radii .\my_radii.csv --density 2500 --seed 42 `
  --output .\results\custom --prefix packing
```

`--buffer` 保留原来的网格间距填充含义，不增加球之间的物理间隙。`--refill-passes` 是累计空遍历上限，成功补填不会消耗这个次数。`--sweeps 0` 可关闭移动，用于诊断。

## 输出

保留四类 CSV 的文件后缀、列名、从 1 开始的 ID、网格词典排序、角点顺序和坐标系平移：

| 文件 | 内容 |
|---|---|
| `*_spheres.csv` | id、x、y、z、radius、diameter、mass |
| `*_summary.csv` | 请求/接受/未放置数量、结束原因、容量标志、总体积、惯性对角线 |
| `*_grid_points.csv` | 有效网格角点的 ID 与坐标 |
| `*_grid_hexahedra.csv` | 八角点连接、球数量、三角面数量 |

命令行保留最终几何、质量、质心和完整惯性张量摘要；新增 `--report` JSON 保存完整 report 字段、惯性张量和各阶段计时。C++ 库通过 `sp::Result` 返回世界坐标状态、输出坐标、质量、总体积、惯性与报告。`capacity_reached` 是保留部分装填结果的正常结束，退出码为 0；非法输入、读写失败返回 1。

## 数值一致性范围

算法和数据格式已保留；不保证跨 MATLAB/C++ 的随机序列或浮点结果逐位相同。相同种子下，C++ 标准库和 MATLAB 的随机数实现不同。`--uniform-tape` 与 `--normal-tape` 支持回放 MATLAB 的原始随机数，文件在计算计时前加载。

规则箱体、容量不足补填及真实 STL 重力压实用例采用严格数值对照。真实佛像长扰动测试中，13,400 次移动顺序和随机数消耗一致；对相同 MATLAB 状态逐次调用 C++ 接触求解，最大接触距离误差为 `7.916e-11`。独立推进两条长扰动轨迹时，最初 `4.44e-16` 的点积舍入差异会被放大，最终坐标最大差约 `0.02916`（半径 `0.625`）。因此长扰动结果不能声称逐坐标相等；基准报告会明确记录差异，并独立检查不重叠、无表面穿透和所有球心位于模型内部。

## 性能与安全设计

- 连续 `std::vector` 存储、稀疏整数网格键、缓存三角面和射线系数。
- DDA 穿越移动路径的全部网格，有限三角形面/边/顶点连续接触，避免穿透薄壁。
- 每次移动立即更新球空间索引；重复三角面使用可复用时间戳去重。
- 保留保守 occupancy 区域分类和边缘处精确奇偶射线测试。
- STL 采用小型严格读取器：支持 ASCII、二进制、以 `solid` 开头的二进制头；校验完整性和有限数值，保持面顺序与退化面，不静默修复输入。
- RAII 管理容器和文件，不使用手动拥有内存的 `new/delete`；非法尺寸/溢出抛出异常。Release 与 AddressSanitizer 测试均已验证。Windows ASAN 检查非法内存访问，不等于形式化的无泄漏证明。
- 当前计算为 CPU 串行，避免改变顺序相关的物理更新；oneAPI 编译器执行优化，没有强行引入 GPU 或共享状态并行。

STL 方案评估了 [stl_reader v2.0](https://github.com/sreiter/stl_reader/tree/v2.0)。其读取流程会删除重复角点构成的退化面，ASCII 解析使用宽松数值转换；这些行为与保留当前输入语义的目标不符，因此采用上述小型严格读取器，不引入运行时网络依赖。

## 复现 MATLAB 对比

本机 MATLAB 使用已安装的 `stlTools/stlRead`，与原程序一致。Python 环境需要 `numpy` 和 `scipy`。

```powershell
# 包含一次预热、三次正式测量、四份 CSV 对照与独立几何校验
.\cpp\scripts\Benchmark.ps1

# 已有 MATLAB 基线时只重跑 C++ 对照
.\cpp\scripts\Benchmark.ps1 -SkipMatlab

# 原有 MATLAB 回归测试：必须显式加入项目绝对路径
matlab -batch "addpath(pwd); addpath(fullfile(pwd,'tests','helpers')); r=runtests('tests'); assertSuccess(r);"

# 严格重放：长扰动的浮点放大将作为失败报告，便于追踪
python .\cpp\benchmarks\compare.py --scenario buddha --strict-replay --verify-only

# 独立检查完整规模结果
python .\cpp\benchmarks\validate_geometry.py `
  --triangles .\cpp\benchmarks\output\buddha\triangles.csv `
  --spheres .\results\Great_Budda_cpp\Great_Budda_packing_spheres.csv
```

计时包括上下文预处理和装填，排除 MATLAB/进程启动、STL 文件导入、CSV 输出、随机数记录及分析器。真实运行的导入与输出时间另存于 JSON。具体实测与局限见 [BENCHMARK_REPORT.md](BENCHMARK_REPORT.md)。
