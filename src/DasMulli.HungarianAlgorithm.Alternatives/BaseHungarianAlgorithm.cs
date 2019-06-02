using System;
using JetBrains.Annotations;
using MathNet.Numerics.LinearAlgebra;

namespace DasMulli
{
    /// <summary>
    /// Class BaseHungarianAlgorithm.
    ///
    /// Implements an assignment algorithm optimizing the global assignment costs.
    /// See https://en.wikipedia.org/wiki/Hungarian_algorithm for an explanation of the algorithm used.
    ///
    /// The implementation is based on @chaowlert's version, Alex Regueiro's implementation at
    /// https://github.com/chaowlert/algorithm/blob/master/Algorithms/HungarianAlgorithm.cs
    /// available under the MIT license
    /// </summary>
    public static class BaseHungarianAlgorithm
    {
        /// <summary>
        /// Finds the assignments with the lowest global assignment cost.
        /// See https://en.wikipedia.org/wiki/Hungarian_algorithm for an explanation of the algorithm used.
        /// </summary>
        /// <param name="assignmentCosts">The assignment costs.</param>
        /// <returns>System.Int32[].</returns>
        /// <exception cref="ArgumentNullException">assignmentCosts</exception>
        /// <exception cref="ArgumentException">This algorithm implementation does not support cost matrices with fewer columns than rows - assignmentCosts</exception>
        public static int[] FindAssignments([NotNull] Matrix<double> assignmentCosts)
        {
            if (assignmentCosts == null)
            {
                throw new ArgumentNullException(nameof(assignmentCosts));
            }

            var costs = assignmentCosts.Clone();

            var rows = costs.RowCount;
            var columns = costs.ColumnCount;

            if (rows > columns)
            {
                throw new ArgumentException("This algorithm implementation does not support cost matrices with fewer columns than rows", nameof(assignmentCosts));
            }

            for (var i = 0; i < rows; i++)
            {
                var min = double.MaxValue;

                for (var j = 0; j < columns; j++)
                {
                    var cost = costs[i, j];

                    if (double.IsNegativeInfinity(cost))
                    {
                        costs[i, j] = cost = double.MinValue;
                    }
                    else if (double.IsPositiveInfinity(cost) || double.IsNaN(cost))
                    {
                        costs[i, j] = cost = double.MaxValue;
                    }

                    min = Math.Min(min, cost);
                }

                if (double.IsInfinity(min))
                {
                    min = double.MinValue;
                }

                for (var j = 0; j < columns; j++)
                {
                    costs[i, j] -= min;
                }
            }

            var masks = new byte[rows, columns];
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

            for (var i = 0; i < rows; i++)
            {
                for (var j = 0; j < columns; j++)
                {
                    if (masks[i, j] == 1)
                    {
                        agentsTasks[i] = j;
                        break;
                    }
                }
            }

            return agentsTasks;

            int RunStep1()
            {
                for (var i = 0; i < rows; i++)
                {
                    for (var j = 0; j < columns; j++)
                    {
                        if (masks[i, j] == 1)
                        {
                            colsCovered[j] = true;
                        }
                    }
                }

                var coveredColsCount = 0;

                for (var j = 0; j < columns; j++)
                {
                    if (colsCovered[j])
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

                    var starColumn = FindStarInRow(masks, columns, zeroLocation.Row);
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
                    var row = FindStarInColumn(masks, rows, path[pathIndex].Column);

                    if (row == -1)
                    {
                        break;
                    }

                    pathIndex++;
                    path[pathIndex] = new Location(row, path[pathIndex - 1].Column);

                    var col = FindPrimeInRow(masks, columns, path[pathIndex].Row);

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

                for (var i = 0; i < rows; i++)
                {
                    for (var j = 0; j < columns; j++)
                    {
                        if (rowsCovered[i])
                        {
                            costs[i, j] += minValue;
                        }

                        if (!colsCovered[j])
                        {
                            costs[i, j] -= minValue;
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
                for (var i = 0; i < rows; i++)
                {
                    for (var j = 0; j < columns; j++)
                    {
                        if (masks[i, j] == 2)
                        {
                            masks[i, j] = 0;
                        }
                    }
                }
            }

            double FindMinimum()
            {
                // Assumes that the costs matrix has been filtered from {-inv,inv,NaN} values

                var minValue = double.MaxValue;

                for (var i = 0; i < rows; i++)
                {
                    for (var j = 0; j < columns; j++)
                    {
                        if (!rowsCovered[i] && !colsCovered[j])
                        {
                            minValue = Math.Min(minValue, costs[i, j]);
                        }
                    }
                }

                return minValue;
            }
        }

        private static int FindPrimeInRow([NotNull] byte[,] masks, int columns, int row)
        {
            if (masks == null)
            {
                throw new ArgumentNullException(nameof(masks));
            }

            for (var j = 0; j < columns; j++)
            {
                if (masks[row, j] == 2)
                {
                    return j;
                }
            }

            return -1;
        }

        private static int FindStarInRow([NotNull] byte[,] masks, int columns, int row)
        {
            if (masks == null)
            {
                throw new ArgumentNullException(nameof(masks));
            }

            for (var j = 0; j < columns; j++)
            {
                if (masks[row, j] == 1)
                {
                    return j;
                }
            }

            return -1;
        }

        private static int FindStarInColumn([NotNull] byte[,] masks, int rows, int column)
        {
            if (masks == null)
            {
                throw new ArgumentNullException(nameof(masks));
            }

            for (var i = 0; i < rows; i++)
            {
                if (masks[i, column] == 1)
                {
                    return i;
                }
            }

            return -1;
        }

        private static bool TryFindZero([NotNull] Matrix<double> costs, [NotNull] bool[] rowsCovered, [NotNull] bool[] colsCovered, out Location zeroLocation)
        {
            if (costs == null)
            {
                throw new ArgumentNullException(nameof(costs));
            }

            if (rowsCovered == null)
            {
                throw new ArgumentNullException(nameof(rowsCovered));
            }

            if (colsCovered == null)
            {
                throw new ArgumentNullException(nameof(colsCovered));
            }

            var storage = costs.AsColumnMajorArray() ?? costs.ToColumnMajorArray();

            if (costs.RowCount >= costs.ColumnCount)
            {
                for (var i = 0; i < costs.RowCount; i++)
                {
                    if (!rowsCovered[i])
                    {
                        for (var j = 0; j < costs.ColumnCount; j++)
                        {
                            if (!colsCovered[j] && storage[j * costs.RowCount + i] <= 0)
                            {
                                zeroLocation = new Location(i, j);
                                return true;
                            }
                        }
                    }
                }
            }
            else
            {
                for (var j = 0; j < costs.ColumnCount; j++)
                {
                    if (!colsCovered[j])
                    {
                        for (var i = 0; i < costs.RowCount; i++)
                        {
                            if (!rowsCovered[i] && storage[j * costs.RowCount + i] <= 0)
                            {
                                zeroLocation = new Location(i, j);
                                return true;
                            }
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
    }
}
