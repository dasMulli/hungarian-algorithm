using System;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using JetBrains.Annotations;
using MathNet.Numerics.LinearAlgebra;

namespace DasMulli
{
    /// <summary>
    /// Class HungarianAlgorithm.
    ///
    /// Implements an assignment algorithm optimizing the global assignment costs.
    /// See https://en.wikipedia.org/wiki/Hungarian_algorithm for an explanation of the algorithm used.
    ///
    /// The implementation is based on @chaowlert's version, Alex Regueiro's implementation at
    /// https://github.com/chaowlert/algorithm/blob/master/Algorithms/HungarianAlgorithm.cs
    /// available under the MIT license
    /// </summary>
    public static class HungarianAlgorithm
    {
        /// <summary>
        /// Finds the assignments with the lowest global assignment cost.
        /// See https://en.wikipedia.org/wiki/Hungarian_algorithm for an explanation of the algorithm used.
        /// </summary>
        /// <param name="assignmentCosts">The assignment costs.</param>
        /// <returns>System.Int32[].</returns>
        /// <exception cref="ArgumentNullException">assignmentCosts</exception>
        /// <exception cref="ArgumentException">This algorithm implementation does not support cost matrices with fewer columns than rows - assignmentCosts</exception>
        public static unsafe int[] FindAssignments([NotNull] Matrix<double> assignmentCosts)
        {
            if (assignmentCosts == null)
            {
                throw new ArgumentNullException(nameof(assignmentCosts));
            }

            var rows = assignmentCosts.RowCount;
            var columns = assignmentCosts.ColumnCount;
            if (rows > columns)
            {
                throw new ArgumentException("This algorithm implementation does not support cost matrices with fewer columns than rows", nameof(assignmentCosts));
            }
            var costs = new FloatStorage(rows, columns);

            for (var i = 0; i < rows; i++)
            {
                var min = float.MaxValue;

                for (var j = 0; j < columns; j++)
                {
                    var cost = (float)assignmentCosts[i, j];

                    if (float.IsNegativeInfinity(cost))
                    {
                        costs[i, j] = cost = float.MinValue;
                    }
                    else if (float.IsPositiveInfinity(cost) || float.IsNaN(cost))
                    {
                        costs[i, j] = cost = float.MaxValue;
                    }
                    else
                    {
                        costs[i, j] = cost;
                    }

                    min = Math.Min(min, cost);
                }

                if (float.IsInfinity(min))
                {
                    min = float.MinValue;
                }

                for (var j = 0; j < columns; j++)
                {
                    costs[i, j] -= min;
                }
            }

            var masks = new ByteStorage(rows, columns);
            var rowsCovered = new bool[rows];
            var colsCovered = new bool[columns];

            for (var i = 0; i < rows; i++)
            {
                for (var j = 0; j < columns; j++)
                {
                    if (costs[i, j] <= 0 && !rowsCovered[i] && !colsCovered[j])
                    {
                        masks[i, j] = 1;
                        rowsCovered[i] = colsCovered[j] = true;
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

            var masksMaxVectorIndex = masks.columnMajorBackingStore.Length -
                                      masks.columnMajorBackingStore.Length % Vector256<byte>.Count;

            if (masksMaxVectorIndex > 0)
            {
                fixed (byte* masksStorePtr = masks.columnMajorBackingStore)
                {
                    var ones = Vector256.Create((byte)1);
                    for (var i = 0; i < masksMaxVectorIndex; i += Vector256<byte>.Count)
                    {
                        var masksVector = Avx.LoadVector256(masksStorePtr + i);
                        var comparisonResult = Avx2.CompareEqual(masksVector, ones);
                        var comparisonMask = (uint)Avx2.MoveMask(comparisonResult);

                        if (comparisonMask == 0)
                        {
                            continue;
                        }

                        var foundIndex = i;
                        while (comparisonMask != 0)
                        {
                            if ((comparisonMask & 1) == 1)
                            {
                                var column = foundIndex / rows;
                                var row = foundIndex % rows;
                                agentsTasks[row] = column;
                            }

                            foundIndex++;
                            comparisonMask >>= 1;
                        }
                    }
                }
            }

            if (masksMaxVectorIndex < masks.columnMajorBackingStore.Length)
            {
                for (var i = masksMaxVectorIndex; i < masks.columnMajorBackingStore.Length; i++)
                {
                    if (masks.columnMajorBackingStore[i] == 1)
                    {
                        var column = i / rows;
                        var row = i % rows;
                        agentsTasks[row] = column;
                    }
                }
            }

            return agentsTasks;

            int RunStep1()
            {
                var maxVectorOffset = rows - rows % Vector256<byte>.Count;

                if (maxVectorOffset > 0)
                {
                    fixed (byte* masksStorePtr = masks.columnMajorBackingStore)
                    {
                        var onesVector = Vector256.Create((byte)1);
                        for (var column = 0; column < columns; column++)
                        {
                            var masksBasePointer = masksStorePtr + column * rows;
                            for (var row = 0; row < maxVectorOffset; row += Vector256<byte>.Count)
                            {
                                var masksRowVector = Avx.LoadVector256(masksBasePointer + row);
                                var comparison = Avx2.CompareEqual(masksRowVector, onesVector);
                                var comparisonMask = Avx2.MoveMask(comparison);

                                if (comparisonMask != 0)
                                {
                                    colsCovered[column] = true;
                                    break;
                                }
                            }
                        }
                    }
                }

                if (maxVectorOffset < rows)
                {
                    for (var column = 0; column < columns; column++)
                    {
                        for (var row = maxVectorOffset; row < rows; row++)
                        {
                            if (masks[row, column] == 1)
                            {
                                colsCovered[column] = true;
                            }
                        }
                    }
                }

                var coveredColsCount = 0;
                var maxColsVectorOffset = columns - columns % Vector256<byte>.Count;

                if (maxColsVectorOffset > 0)
                {
                    var zeroVector = Vector256<byte>.Zero;
                    fixed (bool* colsCoveredPtr = colsCovered)
                    {
                        var colsCoveredBytePtr = (byte*)colsCoveredPtr;
                        for (var column = 0; column < maxColsVectorOffset; column += Vector256<byte>.Count)
                        {
                            var colsCoveredVector = Avx.LoadVector256(colsCoveredBytePtr + column);
                            var zeroMaskVector = Avx2.CompareEqual(colsCoveredVector, zeroVector);
                            var nonZeroMask = ~(uint)Avx2.MoveMask(zeroMaskVector);
                            var count = Popcnt.PopCount(nonZeroMask);
                            coveredColsCount += (int)count;
                        }
                    }
                }

                if (maxColsVectorOffset < columns)
                {
                    for (var column = maxColsVectorOffset; column < columns; column++)
                    {
                        if (colsCovered[column])
                        {
                            coveredColsCount++;
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

                var minValueVector = Vector256.Create(minValue);

                var maxVectorOffset = costs.rowCount - costs.rowCount % Vector256<float>.Count;

                var rowsCoveredFloats = new float[maxVectorOffset];
                for (var row = 0; row < maxVectorOffset; row++)
                {
                    rowsCoveredFloats[row] = rowsCovered[row] ? 1f : 0f;
                }

                fixed (float* rowsCoveredFloatsPtr = rowsCoveredFloats, storagePtr = costs.columnMajorBackingStore)
                    for (var column = 0; column < columns; column++)
                    {
                        var rowBasePtr = storagePtr + column * rows;
                        var columnCoveredVector = Vector256.Create(colsCovered[column] ? 0f : -1f);
                        for (var row = 0; row < maxVectorOffset; row += Vector256<float>.Count)
                        {
                            var rowsCoveredVector = Avx.LoadVector256(rowsCoveredFloatsPtr + row);
                            var multipliers = Avx.Add(rowsCoveredVector, columnCoveredVector);
                            var diff = Avx.Multiply(minValueVector, multipliers);
                            var rowVector = Avx.LoadVector256(rowBasePtr + row);
                            var updatedRow = Avx.Add(rowVector, diff);
                            Avx.Store(rowBasePtr + row, updatedRow);
                        }
                    }

                if (maxVectorOffset < rows)
                {
                    for (var column = 0; column < columns; column++)
                    {
                        for (var row = maxVectorOffset; row < rows; row++)
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
                var storeLength = costs.columnMajorBackingStore.Length;
                var maxVectorOffset = storeLength - storeLength % Vector256<byte>.Count;

                if (maxVectorOffset > 0)
                {
                    var twoVector = Vector256.Create((byte)2);
                    fixed (byte* maskStoragePtr = masks.columnMajorBackingStore)
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
                }

                if (maxVectorOffset < rows)
                {
                    for (var i = maxVectorOffset; i < storeLength; i++)
                    {
                        if (masks.columnMajorBackingStore[i] == 2)
                        {
                            masks.columnMajorBackingStore[i] = 0;
                        }
                    }
                }
            }

            unsafe float FindMinimum()
            {
                // Assumes that the costs matrix has been filtered from {-inv,inv,NaN} values

                var minValue = float.MaxValue;

                var maxVectorOffset = costs.rowCount - costs.rowCount % Vector256<float>.Count;

                var minValueVector = Vector256.Create(minValue);

                var rowsCoveredMaValues = new float[costs.rowCount];
                for (var row = 0; row < costs.rowCount; row++)
                {
                    rowsCoveredMaValues[row] = rowsCovered[row] ? float.MaxValue : 0;
                }

                fixed (float* storagePtr = costs.columnMajorBackingStore)
                fixed (float* rowsCoveredMaValuesPtr = rowsCoveredMaValues)
                    for (var column = 0; column < costs.columnCount; column++)
                    {
                        if (colsCovered[column])
                        {
                            continue;
                        }

                        var rowBasePtr = storagePtr + column * costs.rowCount;

                        for (var row = 0; row < maxVectorOffset; row += Vector256<float>.Count)
                        {
                            var rowVector = Avx.LoadVector256(rowBasePtr + row);
                            var coveredMaxValueVector = Avx.LoadVector256(rowsCoveredMaValuesPtr + row);
                            minValueVector = Avx.Min(Avx.Max(rowVector, coveredMaxValueVector), minValueVector);
                        }
                    }

                var permutedMinValueVector = Avx.Permute2x128(minValueVector, minValueVector, 1);
                var min1 = Avx.Min(minValueVector, permutedMinValueVector);
                var permutedMin1 = Avx.Permute(min1, 0b01_00_11_10);
                var min2 = Avx.Min(min1, permutedMin1);
                var permutedMin2 = Avx.Permute(min2, 0b10_11_00_01);
                var min3 = Avx.Min(min2, permutedMin2);
                minValue = min3.ToScalar();

                if (maxVectorOffset < costs.rowCount)
                {
                    for (var column = 0; column < costs.columnCount; column++)
                    {
                        if (colsCovered[column])
                        {
                            continue;
                        }

                        for (var row = maxVectorOffset; row < costs.rowCount; row++)
                        {
                            if (!rowsCovered[row])
                            {
                                minValue = UnsafeMin(minValue, costs[row, column]);
                            }
                        }
                    }
                }

                return minValue;
            }
        }

        private static int FindPrimeInRow(ByteStorage masks, int row)
        {
            var index = row;
            for (var j = 0; j < masks.columnCount; j++)
            {
                if (masks.columnMajorBackingStore[index] == 2)
                {
                    return j;
                }

                index += masks.rowCount;
            }

            return -1;
        }

        private static int FindStarInRow(ByteStorage masks, int row)
        {
            int index = row;
            for (var j = 0; j < masks.columnCount; j++)
            {
                if (masks.columnMajorBackingStore[index] == 1)
                {
                    return j;
                }

                index += masks.rowCount;
            }

            return -1;
        }

        private static unsafe int FindStarInColumn(ByteStorage masks, int column)
        {
            var rowCount = masks.rowCount;
            var maxVectorOffset = rowCount - rowCount % Vector256<byte>.Count;

            if (maxVectorOffset > 0)
            {
                fixed (byte* storagePtr = masks.columnMajorBackingStore)
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
                            while ((comparisonMask & 1) != 1)
                            {
                                comparisonMask >>= 1;
                                row++;
                            }

                            return row;
                        }
                    }
                }
            }

            if (maxVectorOffset < rowCount)
            {
                for (var row = maxVectorOffset; row < rowCount; row++)
                {
                    if (masks[row, column] == 1)
                    {
                        return row;
                    }
                }
            }

            return -1;
        }

        private static bool TryFindZero(FloatStorage costs, [NotNull] bool[] rowsCovered, [NotNull] bool[] colsCovered, out Location zeroLocation)
        {
            if (colsCovered == null)
            {
                throw new ArgumentNullException(nameof(colsCovered));
            }

            var rowCount = costs.rowCount;
            var columnCount = costs.columnCount;
            var maxVectorOffset = rowCount - rowCount % Vector256<float>.Count;
            var storage = costs.columnMajorBackingStore;
            var zeroVector = Vector256<float>.Zero;

            var coveredMasks = new int[maxVectorOffset / Vector256<float>.Count];
            for (int i = 0; i < maxVectorOffset; i += Vector256<float>.Count)
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

            unsafe
            {
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

                                if ((equality & 0b1) == 0b1)
                                {
                                    zeroLocation = new Location(row, column);
                                    return true;
                                }

                                if ((equality & 0b10) == 0b10)
                                {
                                    zeroLocation = new Location(row + 1, column);
                                    return true;
                                }

                                if ((equality & 0b100) == 0b100)
                                {
                                    zeroLocation = new Location(row + 2, column);
                                    return true;
                                }

                                if ((equality & 0b1000) == 0b1000)
                                {
                                    zeroLocation = new Location(row + 3, column);
                                    return true;
                                }

                                if ((equality & 0b10000) == 0b10000)
                                {
                                    zeroLocation = new Location(row + 4, column);
                                    return true;
                                }

                                if ((equality & 0b100000) == 0b100000)
                                {
                                    zeroLocation = new Location(row + 5, column);
                                    return true;
                                }

                                if ((equality & 0b1000000) == 0b1000000)
                                {
                                    zeroLocation = new Location(row + 6, column);
                                    return true;
                                }

                                if ((equality & 0b10000000) == 0b10000000)
                                {
                                    zeroLocation = new Location(row + 7, column);
                                    return true;
                                }
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

            zeroLocation = new Location(-1, -1);
            return false;
        }

        [SuppressMessage("ReSharper", "InconsistentNaming")]
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

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float UnsafeMin(float left, float right) => left <= right ? left : right;

        private struct FloatStorage
        {
            public readonly int rowCount;
            public readonly int columnCount;
            public readonly float[] columnMajorBackingStore;

            public FloatStorage(int rowCount, int columnCount) : this()
            {
                this.rowCount = rowCount;
                this.columnCount = columnCount;
                columnMajorBackingStore = new float[rowCount * columnCount];
            }

            public float this[int row, int column]
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                get
                {
#if DEBUG
                    Debug.Assert(row >= 0, "Row must be greater than or equal to 0");
                    Debug.Assert(row < rowCount, "Row outside the limit of the current storage");
                    Debug.Assert(column >= 0, "Column must be greater than or equal to 0");
                    Debug.Assert(column < columnCount, "Column outside the limit of the current storage");
#endif
                    return columnMajorBackingStore[column * rowCount + row];
                }
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                set
                {
#if DEBUG
                    Debug.Assert(row >= 0, "Row must be greater than or equal to 0");
                    Debug.Assert(row < rowCount, "Row outside the limit of the current storage");
                    Debug.Assert(column >= 0, "Column must be greater than or equal to 0");
                    Debug.Assert(column < columnCount, "Column outside the limit of the current storage");
#endif
                    columnMajorBackingStore[column * rowCount + row] = value;
                }
            }
        }

        private struct ByteStorage
        {
            public readonly int rowCount;
            public readonly int columnCount;
            public readonly byte[] columnMajorBackingStore;

            public ByteStorage(int rowCount, int columnCount) : this()
            {
                this.rowCount = rowCount;
                this.columnCount = columnCount;
                columnMajorBackingStore = new byte[rowCount * columnCount];
            }

            public byte this[int row, int column]
            {
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                get
                {
#if DEBUG
                    Debug.Assert(row >= 0, "Row must be greater than or equal to 0");
                    Debug.Assert(row < rowCount, "Row outside the limit of the current storage");
                    Debug.Assert(column >= 0, "Column must be greater than or equal to 0");
                    Debug.Assert(column < columnCount, "Column outside the limit of the current storage");
#endif
                    return columnMajorBackingStore[column * rowCount + row];
                }
                [MethodImpl(MethodImplOptions.AggressiveInlining)]
                set
                {
#if DEBUG
                    Debug.Assert(row >= 0, "Row must be greater than or equal to 0");
                    Debug.Assert(row < rowCount, "Row outside the limit of the current storage");
                    Debug.Assert(column >= 0, "Column must be greater than or equal to 0");
                    Debug.Assert(column < columnCount, "Column outside the limit of the current storage");
#endif
                    columnMajorBackingStore[column * rowCount + row] = value;
                }
            }
        }
    }

}
