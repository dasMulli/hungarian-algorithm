using System;
using System.Linq;
using System.Numerics;
using BenchmarkDotNet.Running;
using MathNet.Numerics.LinearAlgebra;

namespace DasMulli.Benchmarks
{
    class Program
    {
        static void Main(string[] args)
        {
            //var rnd = new Random(42);
            //var data = Matrix<double>.Build.Dense(400, 400, (x, y) => rnd.NextDouble() * 100);
            //Console.WriteLine("Starting...");
            //for (int i = 0; i < 3; i++)
            //{
            //    HungarianAlgorithmOptimization4_AvxStep1.FindAssignments(data);
            //    Console.WriteLine("Finished iteration " + i);
            //}
            BenchmarkSwitcher.FromAssembly(typeof(Program).Assembly).Run(args);
        }
    }
}
