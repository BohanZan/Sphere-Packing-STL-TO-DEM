# SpherePacking

SpherePacking 为 DEMBody 生成闭合 STL 几何体内的非重叠球体集合。它接受 STL
文件或内存 mesh，并输出 DEM 球体装配及 CSV 文件。

## 接口与输出

```matlab
[assembly, masses, totalVolume, inertia, report] = ...
    spawnSpheres(model, radii, maxAttempts, buffer, options)
```

- `model` 为闭合 STL 的路径或具有 `vertices`、`faces` 字段的 mesh。
- `assembly` 保持 `4×N`：前三行是球心坐标，第四行是半径；`masses`、
  `totalVolume`、`inertia` 和 `report` 的含义及字段保持不变。
- 指定 `options.outputDirectory` 后，程序保留既有四个 CSV：
  `_spheres.csv`、`_summary.csv`、`_grid_points.csv`、`_grid_hexahedra.csv`。
  球体/点云表头仍为 `id,x,y,z,radius,diameter,mass`。
- `occupancyAcceleration` 在 `spawnSpheres` 中默认启用；它只对保守判定为
  内部/外部的格元提供快捷结果，margin 格元仍使用精确射线奇偶判定。直接以
  四个参数调用 `spBuildContext` 时，该加速层默认关闭。

候选球心通过内部判定后，仍会执行既有的有限三角形距离检查。因此，内部格元的
快捷判定并不表示球体可以穿过 STL 表面。

## 性能验证（2026-09-02）

所有测量使用固定 iron-particle 配置：`inputs/ironParticle/ironParticle.stl`、
100 个半径 `3.0e-6`、`maxAttempts=100`、`buffer=0`、随机种子 53、100 次压缩
sweep、2 次 shake sweep、3 次 refill pass、world 坐标。最终测量关闭 MATLAB
profiler，并写入临时目录。

| 测量 | 三次秒数 | 中位数 | 接受球数 |
| --- | ---: | ---: | ---: |
| 保留的优化前基线 | 126.735856, 128.264100, 146.996799 | 128.264100 | 20, 20, 20 |
| 最终测量 | 72.507578, 73.310399, 75.290296 | 73.310399 | 20, 20, 20 |

最终中位数比保留基线低 42.84%。基线文件未保留当时的 CPU/MATLAB 元数据；因此，
上述比较应理解为相同固定算例的保留基线与本次最终重测，而非跨环境的通用速度承诺。
最终重测环境为 13th Gen Intel(R) Core(TM) i7-13850HX、MATLAB R2025a Update 1
（25.1.0.2973910）。三轮球体体积分数均为 `-0.237247632907`，与保留基线相同；
负号来自该 STL 的有向体积方向。

一次启用/禁用 occupancy 的 profiler 对照保持接受数 20、体积分数完全相同，且将
`spExactPointInside` 调用从 60,745 降至 29,421（减少 51.57%）。此对照在 profiler
开启时运行，故不将其单轮 wall-clock 用作速度比较。

完整 MATLAB 回归为 26/26 通过、0 failed、0 incomplete（11.4865 秒）。四球缩小
profile 记录 19,918 次调用：`verticalTriangleHit` 与 `spPointTriangleDistance` 均为
0 次；`spHashNeighbours`、`spSphereHitsTriangles`、`spExactPointInside` 分别为
391、391、258 次。该缩小 profile 的主要剩余时间包含 CSV 写出和 relaxation，不能
代替完整铁颗粒算例的总耗时测量。

可用下列命令复测最终基准：

```matlab
addpath(pwd); addpath('tests');
profile off; profile clear;
result = benchmarkSpherePacking(3)
```

## 实现范围

优化保留候选顺序、固定种子随机行为、放置、relaxation、refill、CSV 布局和 report
字段。改动仅替换等价的几何查询与 hash/storage 实现：缓存面数据、批量精确射线与
有限三角形判定、保守 occupancy 分流、球邻格收集，以及状态预分配/增量索引。
