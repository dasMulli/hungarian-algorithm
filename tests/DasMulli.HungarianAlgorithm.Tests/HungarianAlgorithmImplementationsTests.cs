using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Threading;
using FluentAssertions;
using MathNet.Numerics.LinearAlgebra;
using Xunit;

namespace DasMulli.Tests
{
    public class HungarianAlgorithmImplementationsTests
    {
        [Theory]
        [MemberData(nameof(ItShallAssignWithMinimumCostTestData))]
        public void ItShallAssignWithMinimumCost(Func<Matrix<double>, int[]> implementation, Matrix<double> assignmentCosts, int[] expectedResult)
        {
            // When
            int[] result = null;
            Exception exception = null;
            var runnerThread = new Thread(() =>
            {
                try
                {
                    result = implementation(assignmentCosts);
                }
                catch (Exception e)
                {
                    exception = e;
                }
            });

            var timeLimit = DateTime.UtcNow + TimeSpan.FromSeconds(2);
            runnerThread.Start();

            bool running = true;
            while (running && timeLimit > DateTime.UtcNow)
            {
                running = !runnerThread.Join(50);
            }

            running.Should().BeFalse("Execution took too long");

            // Then
            exception.Should().BeNull();
            result.Should().NotBeNull();
            result.Should().BeEquivalentTo(expectedResult);
        }

        [Theory]
        [MemberData(nameof(AlgorithmImplementations))]
        public void ItShallThrowExceptionForMatricesWithMoreRowsThanColumns(Func<Matrix<double>, int[]> implementation)
        {
            // Given
            var assignmentCosts =
                Matrix<double>.Build.DenseOfRowArrays(
                    new double[] { 85, 65, 77, 49, 98, 25 },
                    new double[] { 60, 74, 70, 70, 2, 54 },
                    new double[] { 62, 41, 40, 96, 4, 91 },
                    new double[] { 45, 97, 67, 63, 51, 50 },
                    new double[] { 4, 68, 86, 79, 45, 47 },
                    new double[] { 42, 30, 82, 78, 53, 45 },
                    new double[] { 20, 46, 88, 85, 72, 1 },
                    new double[] { 54, 42, 71, 66, 24, 81 });

            // When
            Action when = () => implementation(assignmentCosts);

            // Then
            when.Should().Throw<ArgumentException>().Which.ParamName.Should().StartWith("assignmentCosts");
        }

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
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization4_AvxStep1.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization5_AvxFindMethods.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization6_AvxClearPrimes.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization7_AvxFindMinimum.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization8_AvxStep4.FindAssignments };
                yield return new Func<Matrix<double>, int[]>[] { HungarianAlgorithmOptimization9_AvxAgentStepsResult.FindAssignments };
            }
        }

        public static IEnumerable<object[]> ItShallAssignWithMinimumCostTestData
        {
            get
            {
                foreach (var implementation in AlgorithmImplementations)
                {
                    // r2x2_1
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {14, 15},
                            new double[] {80, 36}),
                        new[] {1, 0}
                    };
                    // r2x2_2
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {78, 80},
                            new double[] {42, 73}),
                        new[] {0, 1}
                    };
                    // r3x3_1
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {80, 7, 93},
                            new double[] {23, 38, 56},
                            new double[] {79, 62, 89}),
                        new[] {2, 1, 0}
                    };
                    // r5X5_1
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {80, 73, 59, 16, 76},
                            new double[] {81, 83, 85, 81, 89},
                            new double[] {85, 30, 57, 10, 41},
                            new double[] {82, 43, 58, 33, 77},
                            new double[] {45, 83, 78, 80, 15}),
                        new[] {2, 0, 1, 4, 3}
                    };
                    // r5X5_2
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {22, 81, 56, 61, 74},
                            new double[] {77, 62, 61, 39, 42},
                            new double[] {90, 21, 28, 11, 30},
                            new double[] {21, 35, 16, 78, 60},
                            new double[] {48, 8, 2, 65, 47}),
                        new[] {2, 0, 1, 3, 4}
                    };
                    // r10x10_1
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {36, 81, 9, 14, 23, 96, 99, 99, 6, 12},
                            new double[] {29, 42, 81, 84, 42, 22, 94, 33, 4, 70},
                            new double[] {30, 15, 50, 83, 26, 20, 79, 67, 7, 20},
                            new double[] {44, 69, 66, 65, 24, 34, 66, 86, 83, 63},
                            new double[] {22, 20, 73, 95, 36, 1, 21, 92, 66, 28},
                            new double[] {56, 60, 98, 45, 33, 3, 83, 51, 52, 38},
                            new double[] {95, 70, 24, 25, 73, 25, 49, 92, 64, 46},
                            new double[] {67, 23, 24, 30, 54, 78, 41, 53, 9, 94},
                            new double[] {79, 16, 2, 79, 94, 43, 28, 65, 50, 45},
                            new double[] {22, 75, 26, 90, 71, 48, 79, 61, 19, 39}),
                        new[] {6, 9, 5, 2, 8, 0, 1, 4, 3, 7}
                    };
                    // r10x10_2
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {58, 27, 76, 49, 23, 81, 75, 73, 60, 14},
                            new double[] {38, 14, 22, 67, 22, 65, 36, 69, 75, 71},
                            new double[] {38, 86, 52, 59, 83, 77, 30, 84, 7, 62},
                            new double[] {5, 28, 15, 25, 65, 44, 49, 3, 35, 54},
                            new double[] {81, 70, 27, 88, 54, 35, 28, 91, 93, 1},
                            new double[] {41, 36, 85, 45, 45, 22, 39, 44, 35, 55},
                            new double[] {32, 63, 96, 75, 94, 4, 43, 68, 46, 63},
                            new double[] {54, 94, 25, 47, 81, 94, 57, 63, 37, 15},
                            new double[] {40, 93, 79, 80, 95, 28, 37, 69, 34, 19},
                            new double[] {25, 47, 12, 83, 53, 16, 61, 37, 45, 48}),
                        new[] {4, 8, 5, 9, 6, 7, 0, 2, 1, 3}
                    };
                    // r3x4_1
                    yield return new[]
                    {
                        implementation[0],
                        Matrix<double>.Build.DenseOfRowArrays(
                            new double[] {71, 78, 13, 45},
                            new double[] {37, 7, 7, 3},
                            new double[] {83, 97, 72, 28}),
                        new[] {2, 1, 3}
                    };
                }
            }
        }
    }
}