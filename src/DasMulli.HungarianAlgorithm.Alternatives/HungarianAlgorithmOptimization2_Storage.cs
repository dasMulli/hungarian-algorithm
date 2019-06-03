using System;
using System.Diagnostics;
using System.Runtime.CompilerServices;
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
    public static class HungarianAlgorithmOptimization2_Storage
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
                for (var row = 0; row < rows; row++)
                {
                    for (var column = 0; column < columns; column++)
                    {
                        if (masks[row, column] == 1)
                        {
                            colsCovered[column] = true;
                        }
                    }
                }

                var coveredColsCount = 0;

                for (var column = 0; column < columns; column++)
                {
                    if (colsCovered[column])
                    {
                        coveredColsCount++;
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
            for (var column = 0; column < masks.ColumnCount; column++)
            {
                if (masks[row, column] == 2)
                {
                    return column;
                }
            }

            return -1;
        }

        private static int FindStarInRow(Storage<byte> masks, int row)
        {
            for (var column = 0; column < masks.ColumnCount; column++)
            {
                if (masks[row, column] == 1)
                {
                    return column;
                }
            }

            return -1;
        }

        private static int FindStarInColumn(Storage<byte> masks, int column)
        {
            for (var row = 0; row < masks.RowCount; row++)
            {
                if (masks[row, column] == 1)
                {
                    return row;
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
#if DEBUG
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
#if DEBUG
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
