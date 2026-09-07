# SpherePacking MATLAB

Static non-overlapping sphere packing in a closed STL domain. The MATLAB
packing flow follows the current sibling `SpherePacking-cpp` implementation.
See [the comparison and verification notes](CPP_ALIGNMENT.md) for the function
mapping, execution order, and floating-point limits.

```matlab
options = struct('initialRelaxationScope','layer', ...
    'maxCompressionSweeps',Inf,'randomSeed',42);
[assembly,masses,volume,inertia,report] = ...
    spawnSpheres(model,radii,1000,0,options);
```

`model` is an STL filename (using the installed `stlRead`) or a struct containing
`vertices` and triangular `faces`. Radii are processed in the supplied order.
`assembly` contains one sphere per column: `[x;y;z;radius]`.

- Omitted/zero buffer selects the maximum input radius. It sets both grid
  padding and the maximum length of each compression or shake move.
- The default initial relaxation scope is `'layer'`; `'batch'` and `'all'`
  select the other C++ policies. Refilling always moves only its new batch.
- The default sweep limit is `Inf`, so energy determines convergence. A finite
  limit raises an error if convergence is not reached. Use `0` to disable motion.
- The default seed is `42`. It repeats MATLAB runs; native C++ random draws
  differ even with the same seed.
- Four CSV files retain the existing sphere, summary, grid-point and
  hexahedron formats, with 17 significant digits for numeric output.

`testRun.m` uses the C++ Buddha preset parameters: 44,000 spheres of radius 0.3.

Run checks from this directory:

```matlab
addpath(pwd,fullfile(pwd,'tests','helpers'));
assertSuccess(runtests(fullfile(pwd,'tests')));
addpath(fullfile(pwd,'tests'));
runCppParity; % Requires the built sibling C++ executable.
```
