#include <iostream>
#include <chrono>
#include <locale.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer.h>
#include <cilk/reducer_opadd.h>
#include <fstream>

const int kPointsNum = 10000000;

const double kIntegrationBorders[2] = { -1, 1 };

/// integrand fucntion
double Func(double x)
{
	return 8 / (2 + 2 * x * x);
}

/// trapezoidal integral calculation
double Integral(int n)
{
	double inegral_sum = 0.0;
	double h = (kIntegrationBorders[1] - kIntegrationBorders[0]) / (n - 1);
	for (auto i = 0; i < n; i++)
		inegral_sum += Func(kIntegrationBorders[0] + i * h) + 
					   Func(kIntegrationBorders[0] + (i + 1)*h);
	return 0.5 * inegral_sum * h;
}

/// parallel trapezoidal integral calculation
double ParallelIntegral(int n)
{
	double h = (kIntegrationBorders[1] - kIntegrationBorders[0]) / (n - 1);
	cilk::reducer_opadd<double> inegral_sum(0.0);
	cilk_for(auto i = 0; i < n; i++)
		inegral_sum += Func(kIntegrationBorders[0] + i * h) + 
					   Func(kIntegrationBorders[0] + (i + 1)*h);
	return 0.5 * inegral_sum->get_value() * h;
}

/// calculation time results logging
void LogTime(char *filename, int n, double time_work)
{
	std::ofstream timework(filename, std::ios::app);
	timework << n << " " << time_work << "\n";
}

int main()
{
	__cilkrts_set_param("nworkers", "4");

	std::chrono::high_resolution_clock::time_point t_begin = std::chrono::high_resolution_clock::now();
	double result = Integral(kPointsNum);
	std::chrono::duration<double> elapsed_time_serial = std::chrono::high_resolution_clock::now() - t_begin;
	std::cout << "Integral = " << result << std::endl;
	std::cout << "Elapsed time = " << elapsed_time_serial.count() << std::endl;

	t_begin = std::chrono::high_resolution_clock::now();
	result = ParallelIntegral(kPointsNum);
	std::chrono::duration<double> elapsed_time_parallel = std::chrono::high_resolution_clock::now() - t_begin;
	std::cout << "Parallel Integral = " << result << std::endl;
	std::cout << "Elapsed time = " << elapsed_time_parallel.count() << std::endl;

	std::cout << "Boost = " << elapsed_time_serial.count() / elapsed_time_parallel.count() << std::endl;

	int points = 1000000000000;
	for (auto i = 100000; i < points; i *= 2)
	{
		t_begin = std::chrono::high_resolution_clock::now();
		result = ParallelIntegral(i);
		elapsed_time_parallel = std::chrono::high_resolution_clock::now() - t_begin;
		LogTime("parallel_result.txt", i, elapsed_time_parallel.count());

	}
	for (auto i = 100000; i < points; i *= 2)
	{
		t_begin = std::chrono::high_resolution_clock::now();
		result = Integral(i);
		elapsed_time_serial = std::chrono::high_resolution_clock::now() - t_begin;
		LogTime("serial_result.txt", i, elapsed_time_serial.count());
	}
	return 0;
}