#include <stdio.h>
#include <ctime>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_vector.h>
#include <cilk/reducer_list.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer.h>

#include <chrono>

using namespace std::chrono;

duration<double> duration_for_s{};
duration<double> duration_for_p{};

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;

/// Функция InitMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void InitMatrix( double** matrix )
{
	for ( int i = 0; i < MATRIX_SIZE; ++i )
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for ( int i = 0; i < MATRIX_SIZE; ++i )
	{
		for ( int j = 0; j <= MATRIX_SIZE; ++j )
		{
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
void SerialGaussMethod( double **matrix, const int rows, double* result )
{
	int k;
	double koef;

	high_resolution_clock::time_point t1, t2;
	t1 = high_resolution_clock::now();
	// прямой ход метода Гаусса
	for ( k = 0; k < rows; ++k )
	{
		//
		for ( int i = k + 1; i < rows; ++i )
		{
			koef = -matrix[i][k] / matrix[k][k];

			for ( int j = k; j <= rows; ++j )
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	t2 = high_resolution_clock::now();
	duration_for_s = (t2 - t1);

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for ( k = rows - 2; k >= 0; --k )
	{
		result[k] = matrix[k][rows];

		//
		for ( int j = k + 1; j < rows; ++j )
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}

	printf("Matrix size:: %dx%d\n", rows, rows);
	printf("Gauss method process time: %lf sec\n", duration_for_s.count());
}


/// Функция ParalGaussMethod() решает СЛАУ методом Гаусса 
/// с учетом параллельной реализации
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
void ParalGaussMethod(double **matrix, const int rows, double* result)
{
	high_resolution_clock::time_point t1, t2;
	t1 = high_resolution_clock::now();
	// прямой ход метода Гаусса
	for (int k = 0; k < rows; ++k)
	{
		//
		cilk_for (int i = k + 1; i < rows; ++i)
		{
			double koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	t2 = high_resolution_clock::now();
	duration_for_p = (t2 - t1);

	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (int k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];

		cilk::reducer_opadd<double> result_k(matrix[k][rows]);
		cilk_for (int j = k + 1; j < rows; ++j)
		{
			//result[k] -= matrix[k][j] * result[j];
			result_k -= matrix[k][j] * result[j];
		}

		//result[k] /= matrix[k][k];
		result[k] = result_k->get_value() / matrix[k][k];
	}

	printf("Matrix size:: %dx%d\n", rows, rows);
	printf("Parallel Gauss method process time: %lf sec\n", duration_for_p.count());
}


int main()
{
	// устанавливаем количество работающих потоков = 4
	__cilkrts_set_param("nworkers", "4");

	srand( (unsigned) time( 0 ) );

	// кол-во строк в матрице, приводимой в качестве примера
	const int test_matrix_lines = 4;

	//double **test_matrix = new double*[test_matrix_lines];
	double **matrix = new double*[MATRIX_SIZE];
	InitMatrix(matrix);

	// цикл по строкам
	/*for (int i = 0; i < test_matrix_lines; ++i )
	{
	//	// (test_matrix_lines + 1)- количество столбцов в тестовой матрице,
	//	// последний столбец матрицы отведен под правые части уравнений, входящих в СЛАУ
		test_matrix[i] = new double[test_matrix_lines + 1];
	}*/

	// массив решений СЛАУ
	//double *result = new double[test_matrix_lines];
	double *result = new double[MATRIX_SIZE];

	// инициализация тестовой матрицы
	/*test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;*/

	//SerialGaussMethod(test_matrix, test_matrix_lines, result );
	SerialGaussMethod(matrix, MATRIX_SIZE, result);
	ParalGaussMethod(matrix, MATRIX_SIZE, result);
	printf("Speed up: %lf\n", duration_for_s.count() / duration_for_p.count());

	/*for (int i = 0; i < test_matrix_lines; ++i )
	{
		delete[] test_matrix[i];
	}*/

	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		delete[] matrix[i];
	}

	/*printf( "Solution:\n" );

	for (int i = 0; i < test_matrix_lines; ++i )
	{
		printf( "x(%d) = %lf\n", i, result[i] );
	}*/

	delete[] matrix;
	delete[] result;

	return 0;
}