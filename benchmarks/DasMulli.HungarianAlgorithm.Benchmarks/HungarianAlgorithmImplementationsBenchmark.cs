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
        public int[] FindAssignmentsOptimization4_AvxStep1() =>
            HungarianAlgorithmOptimization4_AvxStep1.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization5_FindMethodsAvx() =>
            HungarianAlgorithmOptimization5_AvxFindMethods.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization6_AvxClearPrimes() =>
            HungarianAlgorithmOptimization6_AvxClearPrimes.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization7_AvxFindMinimum() =>
            HungarianAlgorithmOptimization7_AvxFindMinimum.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization8_AvxStep4() =>
            HungarianAlgorithmOptimization8_AvxStep4.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignmentsOptimization9_AvxAgentStepsResult() =>
            HungarianAlgorithmOptimization9_AvxAgentStepsResult.FindAssignments(_costs);

        [Benchmark]
        public int[] FindAssignments() => HungarianAlgorithm.FindAssignments(_costs);
    }
}