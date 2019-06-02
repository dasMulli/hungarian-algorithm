using System;
using System.Collections.Generic;
using BenchmarkDotNet.Attributes;
using MathNet.Numerics.LinearAlgebra;

namespace DasMulli.Benchmarks
{
    [RPlotExporter, RankColumn]
    public class HungarianAlgorithmImplementationsBenchmark
    {
        private Matrix<double> _costs;

        [ParamsSource(nameof(CostSizeValues))]
        public int CostSize { get; set; }

        public static IEnumerable<int> CostSizeValues => new[] {50, 100, 250, 500};

        [GlobalSetup]
        public void SetUp()
        {
            var rnd = new Random(42);
            _costs = Matrix<double>.Build.Dense(CostSize, CostSize, (_, __) => rnd.NextDouble() * 100);
        }

        [Benchmark(Baseline = true)]
        public int[] FindAssignmentsBase() => BaseHungarianAlgorithm.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization1_Float() =>
            HungarianAlgorithmOptimization1_Float.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization2_Storage() =>
            HungarianAlgorithmOptimization2_Storage.FindAssignments(_costs);
        [Benchmark]
        public int[] FindAssignmentsOptimization3_AvxFindZero() =>
            HungarianAlgorithmOptimization3_AvxFindZero.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignments() => HungarianAlgorithm.FindAssignments(_costs);
    }
}