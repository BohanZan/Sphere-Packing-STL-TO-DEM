# SpherePacking — Windows C++ / Intel oneAPI

将 MATLAB 实现移植为独立 C++20 库和命令行程序。C++ 源文件位于根目录 `src/`；原 MATLAB 算法和测试集中保存在 `matlab_baseline/`，仅用于基准对照，不参与 C++ 构建。默认构建使用 Intel oneAPI DPC++/C++ Compiler 2025.2，CPU x64、Release 优化和严格浮点运算。

## 项目结构

```text
SpherePacking-cpp/
├── CMakeLists.txt          # 构建配置
├── README.md               # 使用说明
├── include/                # 公共 C++ 头文件（hpp 直接放在此目录）
├── src/                    # C++ 实现及命令行入口
├── tests/                  # C++ 回归测试
│   └── fixtures/           # 测试与参考数据
├── scripts/                # oneAPI 环境、构建、基准运行脚本
├── benchmarks/             # MATLAB/C++ 对照、诊断与性能验证代码
│   ├── trace/              # 仅用于基准的 MATLAB 随机数记录器
│   ├── measured_results/   # 已保存的实测结果快照
│   └── REFERENCES.md       # 算法文献
└── inputs/                 # 输入 STL 模型
```

忽略目录不属于上述源码布局，保持原位：`build/`、`results/`、`tmp/`、`docs/`、`matlab_baseline/`、`Objects/`、`Reference/` 以及已有 `cpp/benchmarks/output/`。`cpp/` 现在仅容纳未迁移的已忽略历史输出，不再包含源码或构建入口。新工作区默认生成 `benchmarks/output/`；当前工作区未创建该目录时，基准脚本会自动复用已有的历史输出，避免搬动数据。

## 编译与运行（PowerShell）

```powershell
Set-Location 'C:\Users\Lenovo\Desktop\Study\AsteroidResearch\DEMBody-move\DEMTools\SpherePacking-cpp'

# 初始化 x64 编译环境、编译并运行回归测试
.\scripts\Build.ps1

# 当前 C++ 预设：请求 44000 个、半径 0.3、重力沿 -Z；迭代至收敛
.\build\release\sphere_packing.exe --preset buddha --report .\results\Great_Budda_cpp\full_report.json
```

构建脚本自动定位 Visual Studio C++ 工具和 oneAPI，优先选择本机 v143/14.44 工具集。运行库会复制到 `.exe` 旁边，运行不需要先进入 oneAPI 命令行。已专门验证普通 PowerShell 以及仅含 Windows 系统路径的环境。

快速运行真实 STL：

```powershell
.\build\release\sphere_packing.exe `
  --model .\inputs\greatBudda\greatBudda.stl `
  --radius 0.625 --count 50 --attempts 60 `
  --shake-sweeps 5 --refill-passes 6 `
  --compression-tolerance 1e-7 --gravity 0 0 -1 --frame world `
  --output .\results\cpp_quick --prefix packing `
  --report .\results\cpp_quick\report.json

# 查看所有选项
.\build\release\sphere_packing.exe --help

# 内存访问检查版本；脚本同时部署 AddressSanitizer DLL
.\scripts\Build.ps1 -Sanitize
.\build\asan\test_packing.exe
```

本机依赖：Intel oneAPI 2025.2、Visual Studio C++ 工具、CMake >= 3.24、Ninja。运行主程序不依赖 MATLAB、Python、网络或第三方 STL 库。MATLAB/NumPy/SciPy 仅用于交叉验证与基准测试。

## 参数与原逻辑

保留有序半径序列、压实与独立侧向扰动、累计空补填遍历上限、任意重力方向、按模型尺度缩放的容差、搜索结束时返回有效部分结果。每次成功生成后立即执行重力压实和扰动压实，再尝试生成下一批；空尝试只累计拒绝次数。

初始装填默认采用候选 A：`--initial-relaxation layer`。进入一个投放层时，程序记录已有球数作为该层的活动起点；该层后续生成成功时，层内累计生成的球共同压实。这个归属按生成批次保留，球沉降到较低位置时仍参与当前层压实。完成整层后才推进活动起点。

可用 `--initial-relaxation all` 选择候选 B，让初始装填阶段所有已有球参与；`--initial-relaxation batch` 保留旧的“仅本次新增球”对照。三者的补填阶段均只移动本轮新补填的球，顶部结束判据也保持相同，以隔离初始活动范围的影响。JSON 报告新增 `initialRelaxationScope`；四类 CSV 格式不变。A 是用户在对照研究中选择的方案，不能宣称其批次划分已被论文唯一确定。

活动范围的选择依据、三种子密实度/耗时对比及完整佛像验证见 [A/B 对照研究](docs/batch-scope-study-2026-09-07.md)。

`--preset buddha` 在 `src/main.cpp` 定义：44000 个球、半径 0.3、重力 -Z、尝试次数 60、每轮扰动 4 次、空补填轮数上限 5、收敛阈值 1e-8、种子 42。输出到 `results/Great_Budda_cpp`；预设后的参数可以覆盖默认值。

不均匀粒径使用无表头的文本/CSV 文件，按顺序给出正半径：

```powershell
.\build\release\sphere_packing.exe --model .\inputs\greatBudda\greatBudda.stl `
  --radii .\my_radii.csv --density 2500 --seed 42 `
  --output .\results\custom --prefix packing
```

`--buffer` 对应论文的 `b_u`：同时控制网格间距 `bL=2*rmax+b_u` 和每次重力/扰动位移上限 `min(contactDistance,b_u)`，不增加物理球间隙。省略或设为 0 表示自动采用输入的最大球半径；正值可显式指定。主程序在构造网格前统一解析该值。

默认仅按势能相对变化收敛，没有固定迭代轮数上限。`--sweeps N` 是显式保护上限；达到上限仍未收敛则退出码为 1，并停止写出本次结果。`--sweeps 0` 是明确关闭移动的诊断模式。`--refill-passes` 是累计空遍历上限，成功补填不会消耗这个次数。

候选球在逆重力方向厚度为 `bL` 的层内生成。论文未明确给出层底的更新公式，本实现从底部开始，每层结束后严格增加 `bL`，保持候选球心高度带相邻连续；空层继续向上推进，以支持不连通几何。球顶高度只用于 Eq. (24) 的顶部结束判断，避免它额外跳过尚未采样的高度带。此锚点规则属于对论文省略细节的实现解释；斜向重力通过主轴直线与斜层相交来采样，不保证斜层内体积均匀分布。球坐标逐个更新，格子列表在每一轮完整的重力或扰动扫描后统一更新。

## 按索引生成反幂律半径

粒径输入仅指定 `--radius 0.3` 时，所有球都采用半径 0.3。球数由 `--count` 指定；省略时使用默认 100，或 `--preset buddha` 的 44000。这种等半径模式与 `--generate-radii inverse-power` 互斥，同时显式指定两者会报错，与参数顺序无关。

反幂律模式的 `sp::generateRadii` 返回含 `count` 项的 `std::vector<double>`，始终覆盖**全部球粒索引 `i=1,...,count`**。先计算严格的反幂律 `raw(i)=C*i^(-p)`，选择 `C` 使**整个截断前数组的中位数**等于指定值，再逐项计算 `radii(i)=clamp(raw(i), minimum, maximum)`。这里的范围只指半径上下界 `--radius-min/max`。偶数个球的中位数是两个中央半径的平均值。生成过程中不抽样、不排序、不打乱索引。

支持两类配置：

- **中位数 + 指数**：指定 `--radius-median` 和 `--radius-exponent`，直接生成；不指定指数或离散程度时，指数默认 1。
- **中位数 + 离散程度**：指定 `--radius-variance` 或 `--radius-stddev`，数值求解指数。指数、方差、标准差三者至多显式指定一个，避免超定约束。方差是总体方差，分母为 N；标准差是其平方根。

中位数和离散程度约束**截断前**的序列，上下界只在最后作用。截断后不会重新拟合，实际中位数和方差可能变化；均值作为计算结果报告，不能再作为第三个独立目标。指数必须非负，0 对应等半径；单个球只能有零方差，两个球的有限反幂律要求方差严格小于中位数的平方。非法或无法数值求解的请求会明确报错。对可表示的截断前数组，还核验实际标准差相对误差约不超过 `1e-8`；例如半径约 1、目标方差 `1e-40` 会被拒绝，因为 double 已无法保存这种微小的半径差异。

下面的 PowerShell 示例生成 200 个半径，覆盖索引 1 到 200，将半径限制在 `[0.2, 0.45]` 后装填真实佛像：

```powershell
.\build\release\sphere_packing.exe --preset buddha `
  --generate-radii inverse-power --count 200 `
  --radius-median 0.3 --radius-variance 0.005 `
  --radius-min 0.2 --radius-max 0.45 `
  --output .\results\variable_radii --prefix packing `
  --report .\results\variable_radii\report.json
```

如需直接控制曲线指数，把 `--radius-variance 0.005` 换成 `--radius-exponent 0.7`；如需标准差入口，可换成 `--radius-stddev 0.07`。`--preset buddha --generate-radii inverse-power` 默认使用 44000 项，预设后的 `--count` 可覆盖数量。球 CSV 的 `radius` 列就是最终半径，正常 `--output` 已包含它，不另写半径文件。

生成器默认参数是中位数 `0.3`、最小半径 `0.1`、最大半径 `0.6`、指数 `1`。这些参数定义在 `include/radii.hpp`；命令行解析位于 `src/main.cpp`。`--radius`、`--radii` 和生成器是不同输入方式，不能把显式 `--radius` 或 `--radii` 与生成器混用。未启用生成器时，原等半径和文件输入方式保持原行为。

作为 C++ 函数调用：

```cpp
#include "radii.hpp"

sp::RadiiOptions input;
input.count = 44000;
input.median = 0.3;
input.variance = 0.01;  // Fit p before clamping; omit an explicit exponent.
input.minimum = 0.15;
input.maximum = 0.6;
sp::RadiiGenerationReport diagnostics;
std::vector<double> radii = sp::generateRadii(input, &diagnostics);
// radii[j] belongs to particle index j+1; all count indices are present.
// Pass this vector directly to sp::pack; build context using max(radii).
```

`diagnostics` 可省略。`sp::radiusStatistics(radii)` 可单独计算均值、中位数、总体方差、标准差、最小值、最大值和数量，不修改原数组。`count` 必须为正且不大于 `2^53-1`，保证转换为 double 后仍可区分每个整数索引。

命令行 `--report` 的 `radiusGeneration` 对象保存数量、实际指数、目标、截断数量，以及 `beforeClamping`、`afterClamping`、`acceptedStatistics` 三组统计。最后一组只统计实际放入模型的球；容量不足时，它们是输入数组的前缀，其统计值可能与完整输入不同。极端指数使截断前半径超出 double 范围时，生成器仍安全输出边界半径，`beforeClamping` 为 `null`；无法表示的方差也用 `null` 标明，报告仍为有效 JSON。零个球被接受时，`acceptedStatistics` 为 `null`。

直接给指数的生成耗时为 O(N)；由方差反求指数需要多次 O(N) 扫描。生成过程不消耗装填随机数，网格和默认 `buffer` 均取**截断后数组**的最大半径。生成耗时单独记为 `radiusGeneration.generationSeconds` 和 `TIMING radii`；原 `computeSeconds`、`totalSeconds` 的阶段口径不变，不包含生成时间。

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

CSV 数据格式保留。2026-09-07 的论文复核发现 MATLAB 基线和旧 C++ 共同存在批次调度、全高度撒球、漏用位移上限等偏差；当前 C++ 已修正，`matlab_baseline/` 保留历史原件。因此不能再把旧 MATLAB 随机数轨迹和计时当作当前算法的等价基线。下面的旧对照结果仅记录修正前版本。

相同种子下，C++ 标准库和 MATLAB 的随机数实现不同。`--uniform-tape` 与 `--normal-tape` 支持外部随机数回放，但必须使用与当前流程一致的记录，旧记录可能耗尽。

规则箱体、容量不足补填及真实 STL 重力压实用例采用严格数值对照。真实佛像长扰动测试中，13,400 次移动顺序和随机数消耗一致；对相同 MATLAB 状态逐次调用 C++ 接触求解，最大接触距离误差为 `7.916e-11`。独立推进两条长扰动轨迹时，最初 `4.44e-16` 的点积舍入差异会被放大，最终坐标最大差约 `0.02916`（半径 `0.625`）。因此长扰动结果不能声称逐坐标相等；基准报告会明确记录差异，并独立检查不重叠、无表面穿透和所有球心位于模型内部。

## 性能与安全设计

- 连续 `std::vector` 存储、稀疏整数网格键、缓存三角面和射线系数。
- DDA 穿越移动路径的全部网格，有限三角形面/边/顶点连续接触，避免穿透薄壁。
- 按 Eq. (22) 限制每次位移，完成整轮扫描后统一更新球空间索引。查询覆盖旧格子列表对应的最大位移，避免漏掉本轮已移动的球；重复三角面使用可复用时间戳去重。
- 势能使用公共尺度归一化后比较相对变化，避免模型单位很大或很小时溢出、下溢导致错误收敛。
- 保留保守 occupancy 区域分类和边缘处精确奇偶射线测试。
- STL 采用小型严格读取器：支持 ASCII、二进制、以 `solid` 开头的二进制头；校验完整性和有限数值，保持面顺序与退化面，不静默修复输入。
- RAII 管理容器和文件，不使用手动拥有内存的 `new/delete`；非法尺寸/溢出抛出异常。Release 与 AddressSanitizer 测试均已验证。Windows ASAN 检查非法内存访问，不等于形式化的无泄漏证明。
- 本机 Intel LLVM 2025.2 的 Windows ASan 存在读取 `catch` 参数崩溃的问题，可用仅抛出/捕获 `std::runtime_error` 的小程序复现，对应 [LLVM PR 159618](https://lists.llvm.org/pipermail/llvm-commits/Week-of-Mon-20250915/1734547.html)。`include/exception.hpp` 仅将异常消息提取函数隔离于插桩之外，算法、几何、I/O 和调用方仍保留检查。CLI 回归测试同时核验退出码和具体错误消息，避免将 ASan 崩溃误判为正常拒绝输入。
- 当前计算为 CPU 串行，避免改变顺序相关的物理更新；oneAPI 编译器执行优化，没有强行引入 GPU 或共享状态并行。

STL 方案评估了 [stl_reader v2.0](https://github.com/sreiter/stl_reader/tree/v2.0)。其读取流程会删除重复角点构成的退化面，ASCII 解析使用宽松数值转换；这些行为与保留当前输入语义的目标不符，因此采用上述小型严格读取器，不引入运行时网络依赖。

## 复现 MATLAB 对比

以下是历史移植版本的对照入口。当前修正后的算法不能直接与旧基线做逐坐标比较；需要先另行建立与论文修正一致的 MATLAB 参考实现。历史加速比不代表当前版本。

本机 MATLAB 使用已安装的 `stlTools/stlRead`，与原程序一致。Python 环境需要 `numpy` 和 `scipy`。

```powershell
# 包含一次预热、三次正式测量、四份 CSV 对照与独立几何校验
.\scripts\Benchmark.ps1

# 已有 MATLAB 基线时只重跑 C++ 对照
.\scripts\Benchmark.ps1 -SkipMatlab

# 原有 MATLAB 回归测试：必须显式加入项目绝对路径
matlab -batch "b=fullfile(pwd,'matlab_baseline'); addpath(b,fullfile(b,'tests','helpers')); r=runtests(fullfile(b,'tests')); assertSuccess(r);"

# 严格重放：长扰动的浮点放大将作为失败报告，便于追踪
python .\benchmarks\compare.py --scenario buddha --strict-replay --verify-only

# 独立检查完整规模结果
$benchmarkOutput = if (Test-Path .\benchmarks\output) { '.\benchmarks\output' } else { '.\cpp\benchmarks\output' }
python .\benchmarks\validate_geometry.py `
  --triangles "$benchmarkOutput\buddha\triangles.csv" `
  --spheres .\results\Great_Budda_cpp\Great_Budda_packing_spheres.csv
```

计时包括上下文预处理和装填，排除 MATLAB/进程启动、STL 文件导入、CSV 输出、随机数记录及分析器。真实运行的导入与输出时间另存于 JSON。历史实测快照见 [benchmarks/measured_results](benchmarks/measured_results/)，算法文献见 [benchmarks/REFERENCES.md](benchmarks/REFERENCES.md)。
