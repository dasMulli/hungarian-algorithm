using System;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using JetBrains.Annotations;
using MathNet.Numerics.LinearAlgebra;

namespace DasMulli
{
    /// <summary>
    /// Class HungarianAlgorithmOptimization2_Storage.
    ///
    /// Implements an assignment algorithm optimizing the global assignment costs.
    /// See https://en.wikipedia.org/wiki/Hungarian_algorithm for an explanation of the algorithm used.
    ///
    /// The implementation is based on @chaowlert's version, Alex Regueiro's implementation at
    /// https://github.com/chaowlert/algorithm/blob/master/Algorithms/HungarianAlgorithm.cs
    /// available under the MIT license
    /// </summary>
    // ReSharper disable once InconsistentNaming
    public static unsafe class HungarianAlgorithmOptimization6_AvxClearPrimes
    {
        /// <summary>
        /// Finds the assignments with the lowest global assignment cost.
        /// See https://en.wikipedia.org/wiki/Hungarian_algorithm for an explanation of the algorithm used.
        /// </summary>
        /// <param name="assignmentCostsDoubles">The assignment costs.</param>
        /// <returns>System.Int32[].</returns>
        /// <exception cref="ArgumentNullException">assignmentCosts</exception>
        /// <exception cref="ArgumentException">This algorithm implementation does not support cost matrices with fewer columns than rows - assignmentCosts</exception>
        public static int[] FindAssignments([NotNull] Matrix<double> assignmentCostsDoubles)
        {
            if (assignmentCostsDoubles == null)
            {
                throw new ArgumentNullException(nameof(assignmentCostsDoubles));
            }

            var rows = assignmentCostsDoubles.RowCount;
            var columns = assignmentCostsDoubles.ColumnCount;

            if (rows > columns)
            {
                throw new ArgumentException("This algorithm implementation does not support cost matrices with fewer columns than rows", nameof(assignmentCostsDoubles));
            }

            var costs = new Storage<float>(rows, columns);

            for (var row = 0; row < rows; row++)
            {
                var min = float.MaxValue;

                for (var column = 0; column < columns; column++)
                {
                    var cost = (float)assignmentCostsDoubles[row, column];

                    if (float.IsNegativeInfinity(cost))
                    {
                        costs[row, column] = cost = float.MinValue;
                    }
                    else if (float.IsPositiveInfinity(cost) || float.IsNaN(cost))
                    {
                        costs[row, column] = cost = float.MaxValue;
                    }
                    else
                    {
                        costs[row, column] = cost;
                    }

                    min = Math.Min(min, cost);
                }

                if (float.IsInfinity(min))
                {
                    min = float.MinValue;
                }

                for (var column = 0; column < columns; column++)
                {
                    costs[row, column] -= min;
                }
            }

            var masks = new Storage<byte>(rows, columns);
            var rowsCovered = new bool[rows];
            var colsCovered = new bool[columns];

            for (var row = 0; row < rows; row++)
            {
                for (var column = 0; column < columns; column++)
                {
                    if (Math.Abs(costs[row, column]) <= 0 && !rowsCovered[row] && !colsCovered[column])
                    {
                        masks[row, column] = 1;
                        rowsCovered[row] = colsCovered[column] = true;
                    }
                }
            }

            ClearCoveredFlags();

            var path = new Location[columns * rows];
            var pathStart = default(Location);
            var step = 1;

            while (step != -1)
            {
                switch (step)
                {
                    case 1:
                        step = RunStep1();
                        break;

                    case 2:
                        step = RunStep2();
                        break;

                    case 3:
                        step = RunStep3();
                        break;

                    case 4:
                        step = RunStep4();
                        break;

                    default:
                        throw new Exception($"Unknown step number {step}");
                }
            }

            var agentsTasks = new int[rows];

            for (var row = 0; row < rows; row++)
            {
                for (var column = 0; column < columns; column++)
                {
                    if (masks[row, column] == 1)
                    {
                        agentsTasks[row] = column;
                        break;
                    }
                }
            }

            return agentsTasks;

            int RunStep1()
            {
                // The covered flags have are reset before this step

                var coveredColsCount = 0;

                if (Avx2.IsSupported && rows >= Vector256<byte>.Count)
                {
                    var maxVectorOffset = rows - rows % Vector256<byte>.Count;

                    if (maxVectorOffset > 0)
                    {
                        var onesVector = Vector256.Create((byte)1);
                        for (var column = 0; column < columns; column++)
                        {
                            var currentColCovered = false;
                            ref var rowRef = ref masks.ColumnMajorBackingStore[column * rows];
                            for (var row = 0; row < maxVectorOffset; row += Vector256<byte>.Count)
                            {
                                var masksRowVector = Unsafe.ReadUnaligned<Vector256<byte>>(ref Unsafe.Add(ref rowRef, row));
                                var comparison = Avx2.CompareEqual(masksRowVector, onesVector);
                                var comparisonMask = Avx2.MoveMask(comparison);

                                if (comparisonMask != 0)
                                {
                                    colsCovered[column] = true;
                                    currentColCovered = true;
                                    coveredColsCount++;
                                    break;
                                }
                            }

                            if (!currentColCovered && maxVectorOffset < rows)
                            {
                                for (var row = maxVectorOffset; row < rows; row++)
                                {
                                    if (Unsafe.Add(ref rowRef, row) == 1)
                                    {
                                        colsCovered[column] = true;
                                        coveredColsCount++;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    for (var column = 0; column < columns; column++)
                    {
                        for (var row = 0; row < rows; row++)
                        {
                            if (masks[row, column] == 1)
                            {
                                colsCovered[column] = true;
                                coveredColsCount++;
                                break;
                            }
                        }
                    }
                }

                if (coveredColsCount == rows)
                {
                    return -1;
                }

                return 2;
            }

            int RunStep2()
            {
                while (true)
                {
                    if (!TryFindZero(costs, rowsCovered, colsCovered, out var zeroLocation))
                    {
                        return 4;
                    }

                    masks[zeroLocation.Row, zeroLocation.Column] = 2;

                    var starColumn = FindStarInRow(masks, zeroLocation.Row);
                    if (starColumn != -1)
                    {
                        rowsCovered[zeroLocation.Row] = true;
                        colsCovered[starColumn] = false;
                    }
                    else
                    {
                        pathStart = zeroLocation;
                        return 3;
                    }
                }
            }

            int RunStep3()
            {
                var pathIndex = 0;
                path[0] = pathStart;

                while (true)
                {
                    var row = FindStarInColumn(masks, path[pathIndex].Column);

                    if (row == -1)
                    {
                        break;
                    }

                    pathIndex++;
                    path[pathIndex] = new Location(row, path[pathIndex - 1].Column);

                    var col = FindPrimeInRow(masks, path[pathIndex].Row);

                    pathIndex++;
                    path[pathIndex] = new Location(path[pathIndex - 1].Row, col);
                }

                ConvertPath(pathIndex + 1);
                ClearCoveredFlags();
                ClearPrimes();

                return 1;
            }

            int RunStep4()
            {
                var minValue = FindMinimum();

                for (var row = 0; row < rows; row++)
                {
                    for (var column = 0; column < columns; column++)
                    {
                        if (rowsCovered[row])
                        {
                            costs[row, column] += minValue;
                        }

                        if (!colsCovered[column])
                        {
                            costs[row, column] -= minValue;
                        }
                    }
                }

                return 2;
            }

            void ClearCoveredFlags()
            {
                Array.Clear(rowsCovered, 0, rowsCovered.Length);
                Array.Clear(colsCovered, 0, colsCovered.Length);
            }

            void ConvertPath(int pathLength)
            {
                for (var i = 0; i < pathLength; i++)
                {
                    var pathItem = path[i];
                    var oldMask = masks[pathItem.Row, pathItem.Column];
                    if (oldMask == 1)
                    {
                        masks[pathItem.Row, pathItem.Column] = 0;
                    }
                    else if (oldMask == 2)
                    {
                        masks[pathItem.Row, pathItem.Column] = 1;
                    }
                }
            }

            void ClearPrimes()
            {
                if (Avx2.IsSupported && rows >= Vector256<byte>.Count)
                {
                    var storeLength = costs.ColumnMajorBackingStore.Length;
                    var maxVectorOffset = storeLength - storeLength % Vector256<byte>.Count;

                    var twoVector = Vector256.Create((byte)2);
                    fixed (byte* maskStoragePtr = masks.ColumnMajorBackingStore)
                    {
                        for (var i = 0; i < maxVectorOffset; i += Vector256<byte>.Count)
                        {
                            var slicePtr = maskStoragePtr + i;
                            var sliceVector = Avx.LoadVector256(slicePtr);
                            var comparison = Avx2.CompareEqual(sliceVector, twoVector);
                            var modifiedRow = Avx2.AndNot(comparison, sliceVector);
                            Avx.Store(slicePtr, modifiedRow);
                        }
                    }

                    if (maxVectorOffset < rows)
                    {
                        for (var i = maxVectorOffset; i < storeLength; i++)
                        {
                            if (masks.ColumnMajorBackingStore[i] == 2)
                            {
                                masks.ColumnMajorBackingStore[i] = 0;
                            }
                        }
                    }
                }
                else
                {
                    for (var row = 0; row < rows; row++)
                    {
                        for (var column = 0; column < columns; column++)
                        {
                            if (masks[row, column] == 2)
                            {
                                masks[row, column] = 0;
                            }
                        }
                    }
                }
            }

            float FindMinimum()
            {
                // Assumes that the costs matrix has been filtered from {-inv,inv,NaN} values

                var minValue = float.MaxValue;

                for (var row = 0; row < rows; row++)
                {
                    for (var column = 0; column < columns; column++)
                    {
                        if (!rowsCovered[row] && !colsCovered[column])
                        {
                            minValue = Math.Min(minValue, costs[row, column]);
                        }
                    }
                }

                return minValue;
            }
        }

        private static int FindPrimeInRow(Storage<byte> masks, int row)
        {
            var index = row;
            for (var j = 0; j < masks.ColumnCount; j++)
            {
                if (masks.ColumnMajorBackingStore[index] == 2)
                {
                    return j;
                }

                index += masks.RowCount;
            }

            return -1;
        }

        private static int FindStarInRow(Storage<byte> masks, int row)
        {
            var index = row;
            for (var j = 0; j < masks.ColumnCount; j++)
            {
                if (masks.ColumnMajorBackingStore[index] == 1)
                {
                    return j;
                }

                index += masks.RowCount;
            }

            return -1;
        }

        private static int FindStarInColumn(Storage<byte> masks, int column)
        {
            if (Avx2.IsSupported && masks.RowCount >= Vector256<byte>.Count)
            {
                var rowCount = masks.RowCount;
                var maxVectorOffset = rowCount - rowCount % Vector256<byte>.Count;

                fixed (byte* storagePtr = masks.ColumnMajorBackingStore)
                {
                    var rowBasePtr = storagePtr + column * rowCount;
                    var ones = Vector256.Create((byte)1);
                    for (var row = 0; row < maxVectorOffset; row += Vector256<byte>.Count)
                    {
                        var masksVector = Avx.LoadVector256(rowBasePtr + row);
                        var comparisonResult = Avx2.CompareEqual(masksVector, ones);
                        var comparisonMask = (uint)Avx2.MoveMask(comparisonResult);

                        if (comparisonMask != 0)
                        {
                            row += (int)Bmi1.TrailingZeroCount(comparisonMask);
                            return row;
                        }
                    }

                    if (maxVectorOffset < rowCount)
                    {
                        rowBasePtr += maxVectorOffset;
                        for (var row = maxVectorOffset; row < rowCount; row++)
                        {
                            if (*rowBasePtr == 1)
                            {
                                return row;
                            }
                            rowBasePtr++;
                        }
                    }
                }
            }
            else
            {
                for (var row = 0; row < masks.RowCount; row++)
                {
                    if (masks[row, column] == 1)
                    {
                        return row;
                    }
                }
            }

            return -1;
        }

        private static bool TryFindZero(Storage<float> costs, [NotNull] bool[] rowsCovered, [NotNull] bool[] colsCovered, out Location zeroLocation)
        {
            if (rowsCovered == null)
            {
                throw new ArgumentNullException(nameof(rowsCovered));
            }

            if (colsCovered == null)
            {
                throw new ArgumentNullException(nameof(colsCovered));
            }

            if (Avx2.IsSupported && costs.RowCount >= Vector256<float>.Count)
            {
                var rowCount = costs.RowCount;
                var columnCount = costs.ColumnCount;
                var storage = costs.ColumnMajorBackingStore;
                var maxVectorOffset = rowCount - rowCount % Vector256<float>.Count;
                var zeroVector = Vector256<float>.Zero;

                var coveredMasks = new int[maxVectorOffset / Vector256<float>.Count];
                for (var i = 0; i < maxVectorOffset; i += Vector256<float>.Count)
                {
                    coveredMasks[i / Vector256<float>.Count] = (rowsCovered[i] ? 0 : 1)
                                      | (rowsCovered[i + 1] ? 0 : 2)
                                      | (rowsCovered[i + 2] ? 0 : 4)
                                      | (rowsCovered[i + 3] ? 0 : 8)
                                      | (rowsCovered[i + 4] ? 0 : 16)
                                      | (rowsCovered[i + 5] ? 0 : 32)
                                      | (rowsCovered[i + 6] ? 0 : 64)
                                      | (rowsCovered[i + 7] ? 0 : 128);
                }

                fixed (float* storagePtr = storage)
                {
                    for (var column = 0; column < columnCount; column++)
                    {
                        if (!colsCovered[column])
                        {
                            var basePtr = storagePtr + rowCount * column;
                            for (int row = 0, rowBatchIndex = 0; row < maxVectorOffset; row += Vector256<float>.Count, rowBatchIndex++)
                            {
                                var rowVector = Avx.LoadVector256(basePtr + row);
                                var comparisonResult = Avx.Compare(rowVector, zeroVector, FloatComparisonMode.OrderedLessThanOrEqualNonSignaling);
                                var equality = (uint)Avx.MoveMask(comparisonResult);

                                if (equality == 0)
                                {
                                    continue;
                                }

                                equality &= (uint)coveredMasks[rowBatchIndex];

                                if (equality == 0)
                                {
                                    continue;
                                }

                                var zeroRow = row + (int)Bmi1.TrailingZeroCount(equality);
                                zeroLocation = new Location(zeroRow, column);
                                return true;
                            }

                            for (var i = maxVectorOffset; i < rowCount; i++)
                            {
                                if (!rowsCovered[i] && storage[column * rowCount + i] <= 0)
                                {
                                    zeroLocation = new Location(i, column);
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                for (var column = 0; column < costs.ColumnCount; column++)
                {
                    if (colsCovered[column]) continue;

                    for (var row = 0; row < costs.RowCount; row++)
                    {
                        if (!rowsCovered[row] && costs.ColumnMajorBackingStore[column * costs.RowCount + row] <= 0)
                        {
                            zeroLocation = new Location(row, column);
                            return true;
                        }
                    }
                }
            }

            zeroLocation = new Location(-1, -1);
            return false;
        }

        private struct Location
        {
            public readonly int Row;
            public readonly int Column;

            public Location(int row, int col)
            {
                Row = row;
                Column = col;
            }
        }

        private struct Storage<T> where T : struct
        {
            public readonly int RowCount;
            public readonly int ColumnCount;
            public readonly T[] ColumnMajorBackingStore;

            public Storage(int rowCount, int columnCount) : this()
            {
                RowCount = rowCount;
                ColumnCount = columnCount;
                ColumnMajorBackingStore = new T[rowCount * columnCount];
            }

            public T this[int row, int column]
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                get
                {
#if DEBUGGER
                    Debug.Assert(row >= 0, "Row must be greater than or equal to 0");
                    Debug.Assert(row < RowCount, "Row outside the limit of the current storage");
                    Debug.Assert(column >= 0, "Column must be greater than or equal to 0");
                    Debug.Assert(column < ColumnCount, "Column outside the limit of the current storage");
#endif
                    return ColumnMajorBackingStore[column * RowCount + row];
                }
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                set
                {
#if DEBUGGER
                    Debug.Assert(row >= 0, "Row must be greater than or equal to 0");
                    Debug.Assert(row < RowCount, "Row outside the limit of the current storage");
                    Debug.Assert(column >= 0, "Column must be greater than or equal to 0");
                    Debug.Assert(column < ColumnCount, "Column outside the limit of the current storage");
#endif
                    ColumnMajorBackingStore[column * RowCount + row] = value;
                }
            }
        }
    }
}
