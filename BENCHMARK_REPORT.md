# Windows / oneAPI 实测报告

测试日期：2026-09-07。源码为本目录 `main-cpp` 工作区，包括用户已有的 MATLAB 未提交修改；这些 MATLAB 文件未被本次移植改写。

## 环境与计时口径

- Windows 11 Pro，10.0.26200。
- Intel Core i7-13850HX，20 核 / 28 逻辑处理器；约 32 GB 内存。
- Intel oneAPI DPC++/C++ Compiler 2025.2.0（Build 20250605），x64、Release、`/fp:strict /Qfma-`。
- MATLAB R2025a Update 1：25.1.0.2973910；使用原有 `stlTools/stlRead` 导入器。
- 每个基准预热一次、正式运行三次，报告中位数。计时包括网格/法线/occupancy 预处理和完整装填阶段；不含 MATLAB/进程启动、STL 导入、CSV 写出、随机数记录或 profiler。C++ 随机回放文件在计时开始前读取。
- CPU 串行 C++；MATLAB 保留现有法线预处理的自动并行策略。桌面机器存在正常调度和热状态波动，原始三次数据全部保留；不将微型用例加速比外推到大规模。

## 时间比较

“回放”列使用与 MATLAB 完全相同的随机数，便于控制工作量；“原生”列是正常 C++ 随机数生成，适合判断实际使用速度。

| 工作负载 | MATLAB 中位数 | C++ 回放中位数 | 回放加速比 | C++ 原生中位数 |
|---|---:|---:|---:|---:|
| 4×4×4 箱体，24 球，r=0.25，压实+扰动 | 0.169904 s | 0.000814 s | 208.8× | 0.000923 s |
| 2×2×2 箱体，请求 20 球，r=0.45，容量不足补填 | 0.025617 s | 0.000187 s | 137.1× | 0.000212 s |
| **greatBudda，89,964 面，50 球，r=0.625，长扰动** | **29.408646 s** | **1.508840 s** | **19.49×** | **1.527190 s（19.26×）** |
| greatBudda，同样 50 球，仅重力压实 | 10.424029 s | 1.415204 s | 7.37× | 1.375003 s |

容量不足用例的 MATLAB/回放均接受 4 球；原生 C++ 因随机序列不同接受 5 球，所以这一行不把原生耗时称为等工作量加速比。其他行均接受目标数量。

长扰动参数：`maxAttempts=60`、`maxCompressionSweeps=200`、`shakeSweeps=5`、`maxRefillPasses=6`、`compressionTolerance=1e-7`、重力 `[0 0 -1]`、seed 42、density 1、world frame、occupancy 开启。纯重力用例改为 sweeps 3、shake 0。它们不是 22,000 球的 MATLAB 全规模计时。

长扰动原始时间：

| 实现 | 第 1 次 | 第 2 次 | 第 3 次 |
|---|---:|---:|---:|
| MATLAB | 29.408646 | 36.455163 | 19.904830 |
| C++ 回放 | 1.531170 | 1.507972 | 1.508840 |
| C++ 原生 | 1.549444 | 1.527190 | 1.482836 |

此前边开发边测量的 C++ 中位数约 1.68 秒；上表使用构建完成后的最终独立重测值。

## 与 MATLAB 的一致性

| 验证 | 结果 |
|---|---|
| 原有 MATLAB 测试 | 80 / 80 通过 |
| C++ Release CTest | 7 / 7 通过（几何、物理、输入输出及 CLI） |
| C++ AddressSanitizer CTest | 7 / 7 通过 |
| 小箱体球数据最大绝对差 | 5.77e-15 |
| 容量不足补填球数据最大绝对差 | 3.77e-15 |
| 真实 STL 纯重力球数据最大绝对差 | 4.97e-14 |
| 四个场景的网格连接、球数、三角面计数 | 与 MATLAB 逐项相同 |
| 相同 MATLAB 状态下的 13,400 次 C++ 接触求解 | 全部误差 < 1e-10；最大 7.916e-11 |
| 四场景独立球间距、有限三角形距离、球心奇偶射线检查 | 通过 |

**长扰动不满足最终坐标数值相等。** 回放同样 417 个均匀、33,000 个正态随机数，移动次序完全相同；第 7 次移动的点积归约产生 `4.44e-16` 差异，随后通过反复接触与侧向移动被放大。最终最大坐标差 `0.0291599`，惯性摘要的最大绝对差约 `0.185415`。这不是可以隐去的误差：测试报告明确标记 `numericalOutputEqualWithin1eMinus8=false`；`--strict-replay` 会为该场景返回失败。

这里保留的是算法、输入/输出接口、物理约束和精度范围，不是 MATLAB 浮点执行过程或随机生成器的逐位仿真。相同原生 seed 也不保证跨语言球坐标或最终装填数量相同。

## 当前 testRun.m 全规模 C++ 运行

`--preset buddha`：请求 22,000 球，r=0.3，gravity `[-1 0 0]`，attempts 40，sweeps 200，shake 5，refill 6，compression tolerance 1e-7，density 1，world frame；C++ seed 42。

| 阶段 | 秒 |
|---|---:|
| STL 导入 | 0.306869 |
| 预处理 | 1.459493 |
| 装填、压实、扰动、补填 | 72.062569 |
| 四份 CSV 输出 | 1.033643 |
| **总计** | **74.862575** |

最终接受 **20,108** 球，未放置 1,892，正常结束原因 `capacity_reached`；完整结果保存在 `results/Great_Budda_cpp/`。所有球经过独立 Python/SciPy 几何校验：最小球间间隙 `-3.109e-15`，最小表面间隙 `-1.610e-15`，模型外球心数量 0。负值处于双精度舍入范围内。

未运行完整 22,000 球 MATLAB baseline，因此不报告该规模的加速倍数。

## 内存与运行库

容器、文件和临时缓冲区均由 RAII 管理；独立代码审查与 ASAN 测试未发现非法内存访问。Windows ASAN 不提供 Linux LSan 那样的完整泄漏检测，这不是形式化的无泄漏证明。

用户截图的缺少 `clang_rt.asan_dynamic-x86_64.dll` 已修复：构建脚本从真实 `compiler/latest` 目录复制非空 DLL，避免 oneAPI 汇总目录中的零长度重定向文件；同时固定 CMake 配置，防止重新识别编译器时丢失 Release/ASAN 开关。已在只保留 Windows 系统 PATH 的进程中启动 ASAN 测试和 Release 主程序，全部通过。

## 复现与原始记录

在项目目录的 PowerShell 中：

```powershell
.\cpp\scripts\Build.ps1
.\cpp\scripts\Benchmark.ps1
.\cpp\scripts\Build.ps1 -Sanitize
```

- `cpp/benchmarks/output/comparison.json`：每次测量、最大误差、接收数量与几何检查。
- `cpp/benchmarks/output/<scenario>/matlab_timing.json`：原始参数、MATLAB 版本与计时。
- `cpp/benchmarks/output/<scenario>/matlab_golden.json`、随机数文件、`matlab/`、`cpp/`：交叉验证记录。
- `tmp/cpp_matlab_tests.log`、`tmp/cpp_matlab_benchmark.log`：MATLAB 运行日志。
- `tmp/trajectory_comparison.json`、`tmp/same_state_summary.txt`：首次舍入分歧与逐状态接触验证。
- `results/Great_Budda_cpp/full_report.json`、`geometry_validation.json`：全规模结果与独立校验。
- `cpp/benchmarks/measured_results/`：本次关键结果的小型快照；`matlab_source_sha256.json` 记录测试用 MATLAB 源文件指纹。

可复现命令和完整参数说明见 [README.md](README.md)。
