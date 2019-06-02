using System.Linq;
using System.Numerics;
using BenchmarkDotNet.Running;

namespace DasMulli.Benchmarks
{
    class Program
    {
        static void Main(string[] args)
        {
            BenchmarkSwitcher.FromAssembly(typeof(Program).Assembly).Run(args);
        }
    }
}
