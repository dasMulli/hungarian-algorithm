using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using FluentAssertions;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace DasMulli.Tests
{
    public class HungarianAlgorithmComparisonTests
    {
        [SuppressMessage("ReSharper", "CoVariantArrayConversion")]
        public static IEnumerable<object[]> AlgorithmImplementations
        {
            get
            {
                yield return new Func<Matrix<double>, int[]>[] { BaseHungarianAlgorithm.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithm.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization1_Float.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization2_Storage.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization3_AvxFindZero.FindAssignments };
            }
        }

        [Theory]
        [MemberData(nameof(AlgorithmImplementations))]
        public void ItShallProcessLargeMatrix(Func<Matrix<double>,int[]> implementation)
        {
            // Given
            var rnd = new Random(42);
            var costs = Matrix<double>.Build.Dense(500, 500, (x, y) => rnd.NextDouble() * 100);

            // When
            implementation(costs);

            // Then no exception is thrown
        }

        [Theory]
        [MemberData(nameof(AlgorithmImplementations))]
        public void ItShallProduceSameResult(Func<Matrix<double>, int[]> implementation)
        {
            // Given
            var rnd = new Random(42);
            var costs = Matrix<double>.Build.Dense(250, 250, (x, y) => rnd.NextDouble() * 100);
            var expected = BaseHungarianAlgorithm.FindAssignments(costs);

            // When
            var actual = implementation(costs);

            // Then no exception is thrown
            expected.Should().BeEquivalentTo(actual);
        }
    }
}
