``` ini

BenchmarkDotNet=v0.11.5, OS=Windows 10.0.17763.529 (1809/October2018Update/Redstone5)
Intel Core i5-8350U CPU 1.70GHz (Kaby Lake R), 1 CPU, 8 logical and 4 physical cores
.NET Core SDK=3.0.100-preview5-011568
  [Host]     : .NET Core 3.0.0-preview5-27626-15 (CoreCLR 4.6.27622.75, CoreFX 4.700.19.22408), 64bit RyuJIT
  DefaultJob : .NET Core 3.0.0-preview5-27626-15 (CoreCLR 4.6.27622.75, CoreFX 4.700.19.22408), 64bit RyuJIT


```
|                   Method | CostSize |     Mean |     Error |    StdDev | Ratio | RatioSD |
|------------------------- |--------- |---------:|----------:|----------:|------:|--------:|
|        **BaseFindStarInRow** |       **50** | **3.076 ns** | **0.0963 ns** | **0.0989 ns** |  **1.00** |    **0.00** |
| FindStarInRowDirectIndex |       50 | 3.062 ns | 0.1647 ns | 0.1961 ns |  1.01 |    0.09 |
|                          |          |          |           |           |       |         |
|        **BaseFindStarInRow** |      **100** | **5.775 ns** | **0.1472 ns** | **0.1965 ns** |  **1.00** |    **0.00** |
| FindStarInRowDirectIndex |      100 | 5.266 ns | 0.0983 ns | 0.0872 ns |  0.91 |    0.04 |
|                          |          |          |           |           |       |         |
|        **BaseFindStarInRow** |      **250** | **6.328 ns** | **0.0543 ns** | **0.0453 ns** |  **1.00** |    **0.00** |
| FindStarInRowDirectIndex |      250 | 6.018 ns | 0.1633 ns | 0.1528 ns |  0.95 |    0.02 |
|                          |          |          |           |           |       |         |
|        **BaseFindStarInRow** |      **500** | **5.331 ns** | **0.1862 ns** | **0.1650 ns** |  **1.00** |    **0.00** |
| FindStarInRowDirectIndex |      500 | 4.394 ns | 0.0876 ns | 0.0819 ns |  0.82 |    0.03 |
