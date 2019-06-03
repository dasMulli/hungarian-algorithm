using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using BenchmarkDotNet.Attributes;

namespace DasMulli.Benchmarks
{
    public class FindStarInRowBenchmark
    {
        private Storage<byte> _masks;

        [ParamsSource(nameof(CostSizeValues))]
        public int CostSize { get; set; }

        public static IEnumerable<int> CostSizeValues => new[] { 50, 100, 250, 500 };


        [GlobalSetup]
        public void SetUp()
        {
            var rnd = new Random(42);
            _masks = new Storage<byte>(CostSize, CostSize);
            for (var i = 0; i < CostSize * CostSize; i++)
            {
                _masks.ColumnMajorBackingStore[i] = rnd.NextDouble() switch
                {
                    var rndVal when rndVal > 0.9 => (byte)2,
                    var rndVal when rndVal > 0.7 => (byte)1,
                    _ => (byte)0
                };
            }
        }

        [Benchmark(Baseline = true)]
        public int BaseFindStarInRow()

        {
            var row = 42;

            for (var column = 0; column < _masks.ColumnCount; column++)
            {
                if (_masks[row, column] == 1)
                {
                    return column;
                }
            }

            return -1;
        }

        [Benchmark]
        public int FindStarInRowDirectIndex()

        {
            var row = 42;

            var index = row;
            for (var j = 0; j < _masks.ColumnCount; j++)
            {
                if (_masks.ColumnMajorBackingStore[index] == 1)
                {
                    return j;
                }

                index += _masks.RowCount;
            }

            return -1;
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
