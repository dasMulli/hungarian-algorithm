using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using FluentAssertions;
using MathNet.Numerics.LinearAlgebra;
using Xunit;
using Xunit.Abstractions;

using FindAssignmentsImpl = System.Func<MathNet.Numerics.LinearAlgebra.Matrix<double>, int[]>;

namespace DasMulli.Tests
{
    public class HungarianAlgorithmComparisonTests
    {
        private readonly ITestOutputHelper output;

        public HungarianAlgorithmComparisonTests(ITestOutputHelper output)
        {
            this.output = output;
        }

        [SuppressMessage("ReSharper", "CoVariantArrayConversion")]
        public static IEnumerable<object[]> AlgorithmImplementations
        {
            get
            {
                yield return new object[] { (FindAssignmentsImpl)BaseHungarianAlgorithm.FindAssignments, nameof(BaseHungarianAlgorithm) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithm.FindAssignments, nameof(HungarianAlgorithm) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization1_Float.FindAssignments, nameof(HungarianAlgorithmOptimization1_Float) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization2_Storage.FindAssignments, nameof(HungarianAlgorithmOptimization2_Storage) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization3_AvxFindZero.FindAssignments, nameof(HungarianAlgorithmOptimization3_AvxFindZero) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization4_AvxStep1.FindAssignments, nameof(HungarianAlgorithmOptimization4_AvxStep1) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization5_AvxFindMethods.FindAssignments, nameof(HungarianAlgorithmOptimization5_AvxFindMethods) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization6_AvxClearPrimes.FindAssignments, nameof(HungarianAlgorithmOptimization6_AvxClearPrimes) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization7_AvxFindMinimum.FindAssignments, nameof(HungarianAlgorithmOptimization7_AvxFindMinimum) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization8_AvxStep4.FindAssignments, nameof(HungarianAlgorithmOptimization8_AvxStep4) };
                yield return new object[] { (FindAssignmentsImpl)HungarianAlgorithmOptimization9_AvxAgentStepsResult.FindAssignments, nameof(HungarianAlgorithmOptimization9_AvxAgentStepsResult) };
            }
        }

        [Theory]
        [MemberData(nameof(AlgorithmImplementations))]
        public void ItShallProcessLargeMatrix(Func<Matrix<double>,int[]> implementation, string name)
        {
            // Given
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var rnd = new Random(42);
            var costs = Matrix<double>.Build.Dense(500, 500, (x, y) => rnd.NextDouble() * 100);
            stopWatch.Stop();
            output.WriteLine("Constructing matrix took {0} ms", stopWatch.ElapsedMilliseconds);

            // When
            stopWatch.Restart();
            implementation(costs);
            stopWatch.Stop();
            output.WriteLine("Running implementation took {0} ms", stopWatch.ElapsedMilliseconds);

            // Then no exception is thrown
        }

        [Theory]
        [MemberData(nameof(AlgorithmImplementations))]
        public void ItShallProduceSameResult(Func<Matrix<double>, int[]> implementation, string name)
        {
            // Given
            var rnd = new Random(42);
            var costs = Matrix<double>.Build.Dense(250, 250, (x, y) => rnd.NextDouble() * 100);
            var expected = BaseHungarianAlgorithm.FindAssignments(costs.Clone());

            // When
            var actual = implementation(costs.Clone());

            // Then no exception is thrown
            expected.Should().BeEquivalentTo(actual);
        }
    }
}
