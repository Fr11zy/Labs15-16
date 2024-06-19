#include <iostream>
#include <string>
#include "Matrix.h"

int main() {   
    //Для файлика
    std::string inputFile;
    std::cout << "Matrix input file > ";
    std::cin >> inputFile;
    std::ifstream inmatrix(inputFile);

    //Записываем матрицы с файлика и консоли
    Matrix<double> matr_from_file;
    inmatrix >> matr_from_file;

    Matrix<double> matr_from_console;
    std::cin >> matr_from_console;

    //Все записыывается в файлик
    std::string outputFile;
    std::cout << "Matrix output file > ";
    std::cin >> outputFile;
    std::ofstream outresult(outputFile);

    Matrix matr1=Matrix(matr_from_console);
    matr1.swap_N(0,1);

    Matrix matr2=Matrix(matr_from_file);
    matr2.multi_N(0,5);

    Matrix matr3=Matrix(matr_from_console);
    matr3.multi_add(0,1,3);


    //Вывод
    std::cout << "Изначальная матрица" << std::endl;
    std::cout << matr_from_file << std::endl;
    outresult << matr_from_console << std::endl;
    
    std::cout << "Сложение двух матриц" << std::endl;
    std::cout << matr1+matr2 << std::endl;
    outresult << matr1+matr2 << std::endl;

    std::cout << "Элементарное преобразование 1 типа" << std::endl;
    std::cout << matr1 << std::endl;
    outresult << matr1 << std::endl;

    std::cout << "Элементарное преобразование 2 типа" << std::endl;
    std::cout << matr2 << std::endl;
    outresult << matr2 << std::endl;

    std::cout << "Элементарное преобразование 3 типа" << std::endl;
    std::cout << matr3 << std::endl;
    outresult << matr3 << std::endl;

    std::cout << "Обратная матрица" << std::endl;
    std::cout << !matr_from_file << std::endl;
    outresult << !matr_from_file << std::endl;

    Matrix<double> sum_result_parallel = matr_from_file.sum_parallel(matr_from_console);
    Matrix<double> subtract_result_parallel = matr_from_file.subtract_parallel(matr_from_console);
    Matrix<double> multiply_result_parallel = matr_from_file.multiply_parallel(matr_from_console);

    //Параллельное сложение с использованием std::async
    auto future_F = matr_from_file.sum_async(matr_from_console,2);
    Matrix<double> matrF = future_F.get();
    std::cout << "Сложение двух матриц(async)" << std::endl;
    std::cout << matrF << std::endl;
    outresult << matrF << std::endl; 

    std::cout << "Сложение двух матриц (параллельно)" << std::endl;
    std::cout << sum_result_parallel << std::endl;
    outresult << sum_result_parallel << std::endl;

    return 0;
}
