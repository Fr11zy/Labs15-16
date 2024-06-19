#ifndef MATRIX_H_
#define MATRIX_H_

#include <fstream>
#include <thread>
#include <future>
#include <chrono>
#include <stdexcept>
#include <vector>

#define Pogr 1e-9

template<typename T>
class Matrix {
    private:
        int N;
        int M;
        T **A;

    public:
        //Конструкторы и деструктор
        ~Matrix()
        {
            for (int indexN=0;indexN<N;indexN++)
            {
                delete [] A[indexN];
            }
            delete [] A;
        }

        Matrix() : N(0), M(0) {}

        Matrix(int newN, int newM)
        {
            if ((newN<=0) || (newM<=0))
                throw std::invalid_argument("Число строк и столбцов должно быть положительным числом");
            this->N=newN;
            this->M=newM;
            A=new T *[N];
            for (int indexN=0;indexN<N;indexN++)
            {
                this->A[indexN]=new T [M];
                for (int indexM=0;indexM<M;indexM++)
                {
                    A[indexN][indexM]=0;
                }
            }
        }

        Matrix(int newN, int newM, T **newA)
        {
            if ((newN<=0) || (newM<=0))
                throw std::invalid_argument("Число строк и столбцов должно быть положительным числом");
            this->N=newN;
            this->M=newM;
            A=new T *[N];
            for (int indexN=0;indexN<N;indexN++)
            {
                for (int indexM=0;indexM<M;indexM++)
                {
                    A[indexN][indexM]=std::abs(newA[indexN][indexM])<Pogr ? 0 : newA[indexN][indexM];
                }
            }    
        }

        Matrix(const Matrix& other)
        {
            this->N=other.N;
            this->M=other.M;
            A=new T *[N];
            for (int indexN=0;indexN<N;indexN++)
            {
                this->A[indexN]=new T [M];
                for (int indexM=0;indexM<M;indexM++)
                {
                    A[indexN][indexM]=other.A[indexN][indexM];
                }
            }
        }
        explicit Matrix(int N) : Matrix(N,N) {}// для того чтобы создать квадратную матрицу
        

        static Matrix matr_one(int N)// чтобы создать единичную матрицу
        {
            Matrix result(N);
            for (int i=0;i<N;i++)
            {
                result.set(i,i,1);
            }
            return result;
        }

        void changeSize(int newN,int newM)// чтобы изменять размер матрицы
        {
            for (int indexN=0;indexN<N;indexN++)
                delete [] A[indexN];
            delete [] A;
            this->N=newN;
            this->M=newM;
            A=new T* [N];
            for (int indexN=0;indexN<N;indexN++)
            {
                A[indexN]=new T [M];
                for (int indexM=0;indexM<M;indexM++)
                {
                    A[indexN][indexM]=0;
                }
            }

        }

        //Логические операторы
        bool operator==(const Matrix& other) const
        {
            if ((this->getN()!=other.getN()) || (this->getM()!=other.getM())) 
                return false;
            for (int indexN=0;indexN<this->getN();indexN++)
            {
                for (int indexM=0;indexM<this->getM();indexM++)
                {
                    if ((std::abs((this->get(indexN,indexM))-(other.get(indexN,indexM))))>Pogr)
                        return false;
                }
            }
            return true;    
        }
        bool operator!=(const Matrix& other) const
        {
            return !(*this==other);
        }
        bool operator==(const double &lambda) const
        {
            return *this==Matrix::matr_one(this->getN())*lambda;
        }
        bool operator!=(const double &lambda) const
        {
            return !(*this==lambda);
        }

        //Чтобы принимать значения и устанавливать значения
        const T &get(int N, int M) const
        {
            if (!(this->isIndexesOK(N,M)))
                throw std::invalid_argument("Максимальное количество строк - "+std::to_string(this->getN()-1)+","+"Максимальное количество столбцов - "+std::to_string(this->getM()-1));
            return this->A[N][M];    
        }
        void set(int N, int M, T newelem)
        {
            if (!(this->isIndexesOK(N,M)))
                throw std::invalid_argument("Максимальное количество строк - "+std::to_string(this->getN()-1)+","+"Максимальное количество столбцов - "+std::to_string(this->getM()-1));
            this->A[N][M]=std::abs(newelem)< Pogr? 0: newelem; 
        }
        int getN() const 
        {
            return N;
        }
        int getM() const 
        {
            return M;
        }
        bool isIndexesOK(unsigned long long N, unsigned long long M) const 
        {
            return N < this->getN() &&
                M < this->getM();
        }

        //Арифметические операции
        Matrix operator+(const Matrix& other) const
        {
            if (((this->getN())!=(other.getN())) || ((this->getM())!=(other.getM())))
                throw std::invalid_argument("Размеры матриц должны быть одинаковыми");
            Matrix result(this->getN(),this->getM());
            for (int indexN=0;indexN<this->getN();indexN++)
            {
                for (int indexM=0;indexM<this->getM();indexM++)
                {
                    result.set(indexN,indexM,((this->get(indexN,indexM))+(other.get(indexN,indexM))));
                }
            }
            return result;
        }

        Matrix operator-(const Matrix& other) const
        {
            if (((this->getN())!=(other.getN())) || ((this->getM())!=(other.getM())))
                throw std::invalid_argument("Размеры матриц должны быть одинаковыми");
            Matrix result(this->getN(),this->getM());
            for (int indexN=0;indexN<this->getN();indexN++)
            {
                for (int indexM=0;indexM<this->getM();indexM++)
                {
                    result.set(indexN,indexM,((this->get(indexN,indexM))-(other.get(indexN,indexM))));
                }
            }
            return result;    
        }

        Matrix operator*(const Matrix& other) const
        {
            if ((this->getN())!=(this->getM()))
                throw std::invalid_argument("Количество столбцов 1-ой матрицы должно совпадать c количеством строк 2-ой матрицы");
            Matrix result(this->getN(),other.getM());
            for (int indexN=0;indexN<N;indexN++)
            {
                for (int indexM=0;indexM<other.M;indexM++)
                {
                    T sigma_for_elem=0;
                    for (int i=0;i<M;i++)
                    {
                        sigma_for_elem+=(this->get(indexN,i))*(other.get(i,indexM));
                    }
                    result.set(indexN,indexM,sigma_for_elem);
                }
            }
            return result;
        }

        Matrix operator*(const double &lambda) const
        {
            Matrix result(this->getN(),this->getM());
            for (int indexN=0;indexN<this->getN();indexN++)
            {
                for (int indexM=0;indexM<this->getM();indexM++)
                {
                    result.set(indexN,indexM,(this->get(N,M))*lambda);
                }
            }
            return result;
        }

        //Элементарные проеобразования
        void swap_N(int N1,int N2)
        {
            if (!(this->isIndexesOK(N1,0)) || !(this->isIndexesOK(N2,0)))
                throw std::invalid_argument("Максимальное число строк"+std::to_string(this->getN()-1));
            std::swap(this->A[N1],this->A[N2]);
        }
        void multi_N(int N,double lambda)
        {   
            if (!(this->isIndexesOK(N,0)))
                throw std::invalid_argument("Максимальное число строк"+std::to_string(this->getN()-1));
            if (lambda==0)
                throw std::invalid_argument("Число lambda должно быть не равно 0");
            for (int indexM=0;indexM<this->getM();indexM++)
            {
                this->set(N,indexM,(this->get(N,indexM))*lambda);
            }
        }
        void multi_add(int N1,int N2, double lambda)//N1*lambda+N2
        {
            if (!(this->isIndexesOK(N1,0)) || !(this->isIndexesOK(N2,0)))
                throw std::invalid_argument("Максимальное число строк"+std::to_string(this->getN()-1));
            if (lambda==0)
                throw std::invalid_argument("Число lambda должно быть не равно 0");
            if (N1==N2)
                throw std::invalid_argument("Номера строк должен быть различными");
            for (int indexM=0;indexM<this->getM();indexM++)
            {
                this->set(N2,indexM,(this->get(N1,indexM))*lambda+(this->get(N2,indexM)));
            }
        }

        T determinant() const {
            if (this->N != this->M)
                throw std::invalid_argument("Матрица должна быть квадратной");

            Matrix<T> lu(*this);
            T det = 1;
            for (int k=0; k<N; k++) {
                int pivot = k;
                for (int i=k+1; i<N; i++) {
                    if (std::abs(lu.get(i,k)) > std::abs(lu.get(pivot,k))) {
                        pivot = i;
                    }  
                }
                if (std::abs(lu.get(pivot,k)) < Pogr) {
                    return 0;
                }
                if (k != pivot) {
                    lu.swap_N(k, pivot);
                    det=-det;
                }
                det *=lu.get(k,k);
                for (int i=k+1; i<this->getN(); i++) {
                    lu.set(i,k,lu.get(i,k)/lu.get(k,k));
                    for (int j = k+1; j<this->getN(); j++) {
                        lu.set(i,j,lu.get(i,j)-lu.get(i,k)*lu.get(k,j));
                    }
                }
            }
            return det;
        }

        Matrix operator!() const {
            if (this->getN() != this->getM()) 
                throw std::invalid_argument("Матрица должна быть квадратной");

            Matrix<T> aug(N, 2*N);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                    aug.set(i,j,this->get(i,j));
                aug.set(i, N+1, 1);
            }

            for (int i = 0; i < N; i++)
            {
                int pivot = i;
                for (int j = i + 1; j < N; j++)
                {
                    if (std::abs(aug.A[j][i]) > std::abs(aug.A[pivot][i]))
                        pivot = j;
                }
                if (std::abs(aug.A[pivot][i]) < Pogr)
                    throw std::invalid_argument("Матрица вырождена и не имеет обратной");
                if (i != pivot)
                    aug.swap_N(i, pivot);

                T pivotValue = aug.A[i][i];
                for (int j = 0; j < 2 * N; j++)
                    aug.A[i][j] /= pivotValue;

                for (int j = 0; j < N; j++)
                {
                    if (i != j)
                    {
                        T factor = aug.A[j][i];
                        for (int k = 0; k < 2 * N; k++)
                            aug.A[j][k] -= factor * aug.A[i][k];
                    }
                }
            }

            Matrix<T> inv(N, N);
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    inv.A[i][j] = aug.A[i][N + j];
            return inv;
        }

        // Потоки task 1
        Matrix sum_parallel(const Matrix& other) const
        {
            if (N != other.N || M != other.M)
            {
                throw std::invalid_argument("Размеры матриц должны быть одинаковыми");
            }

            Matrix<T> result(N, M);
            int num_threads = std::thread::hardware_concurrency();
            int rows_per_thread = N / num_threads;

            auto add_blocks = [&](int start_row, int end_row) {
                for (int i = start_row; i < end_row; ++i)
                {
                    for (int j = 0; j < M; ++j)
                    {
                        result.set(i, j, this->get(i, j) + other.get(i, j));
                    }
                }
            };

            std::vector<std::thread> threads;
            for (int i = 0; i < num_threads; ++i)
            {
                int start_row = i * rows_per_thread;
                int end_row = (i == num_threads - 1) ? N : (i + 1) * rows_per_thread;
                threads.push_back(std::thread(add_blocks, start_row, end_row));
            }

            for (auto& t : threads)
            {
                t.join();
            }

            return result;
        }

        Matrix subtract_parallel(const Matrix& other) const
        {
            if (N != other.N || M != other.M)
            {
                throw std::invalid_argument("Размеры матриц должны быть одинаковыми");
            }

            Matrix<T> result(N, M);
            int num_threads = std::thread::hardware_concurrency();
            int rows_per_thread = N / num_threads;

            auto subtract_blocks = [&](int start_row, int end_row) {
                for (int i = start_row; i < end_row; ++i)
                {
                    for (int j = 0; j < M; ++j)
                    {
                        result.set(i, j, this->get(i, j) - other.get(i, j));
                    }
                }
            };

            std::vector<std::thread> threads;
            for (int i = 0; i < num_threads; ++i)
            {
                int start_row = i * rows_per_thread;
                int end_row = (i == num_threads - 1) ? N : (i + 1) * rows_per_thread;
                threads.push_back(std::thread(subtract_blocks, start_row, end_row));
            }

            for (auto& t : threads)
            {
                t.join();
            }

            return result;
        }

        Matrix multiply_parallel(const Matrix& other) const
        {
            if (M != other.N)
            {
                throw std::invalid_argument("Количество столбцов первой матрицы должно совпадать с количеством строк второй матрицы");
            }

            Matrix<T> result(N, other.M);
            int num_threads = std::thread::hardware_concurrency();
            int rows_per_thread = N / num_threads;

            auto multiply_blocks = [&](int start_row, int end_row) {
                for (int i = start_row; i < end_row; ++i)
                {
                    for (int j = 0; j < other.M; ++j)
                    {
                        T sum = 0;
                        for (int k = 0; k < M; ++k)
                        {
                            sum += this->get(i, k) * other.get(k, j);
                        }
                        result.set(i, j, sum);
                    }
                }
            };

            std::vector<std::thread> threads;
            for (int i = 0; i < num_threads; ++i)
            {
                int start_row = i * rows_per_thread;
                int end_row = (i == num_threads - 1) ? N : (i + 1) * rows_per_thread;
                threads.push_back(std::thread(multiply_blocks, start_row, end_row));
            }

            for (auto& t : threads)
            {
                t.join();
            }

            return result;
        }

        // Потоки task 2
        std::future<Matrix<T>> sum_async(const Matrix<T>& other, unsigned int blocks) const
        {
            if (N != other.N || M != other.M)
            {
                throw std::invalid_argument("Размеры матриц должны быть одинаковыми");
            }

            Matrix<T> result(N, M);
            auto add_blocks = [&](int start_row, int end_row) {
                for (int i = start_row; i < end_row; ++i)
                {
                    for (int j = 0; j < M; ++j)
                    {
                        result.set(i, j, this->get(i, j) + other.get(i, j));
                    }
                }
            };

            std::vector<std::future<void>> futures;
            for (int start_row = 0; start_row < N; start_row += static_cast<int>(blocks))
            {
                int end_row = std::min(start_row + static_cast<int>(blocks), N);
                futures.push_back(std::async(std::launch::async, add_blocks, start_row, end_row));
            }

            for (auto& f : futures)
            {
                f.get();
            }

            return std::async(std::launch::deferred, [result]() { return result; });
        }

        std::future<Matrix<T>> subtract_async(const Matrix<T>& other, unsigned int blocks) const
        {
            if (N != other.N || M != other.M)
            {
                throw std::invalid_argument("Размеры матриц должны быть одинаковыми");
            }

            Matrix<T> result(N, M);
            auto subtract_blocks = [&](int start_row, int end_row) {
                for (int i = start_row; i < end_row; ++i)
                {
                    for (int j = 0; j < M; ++j)
                    {
                        result.set(i, j, this->get(i, j) - other.get(i, j));
                    }
                }
            };

            std::vector<std::future<void>> futures;
            for (int start_row = 0; start_row < N; start_row += static_cast<int>(blocks))
            {
                int end_row = std::min(start_row + static_cast<int>(blocks), N);
                futures.push_back(std::async(std::launch::async, subtract_blocks, start_row, end_row));
            }

            for (auto& f : futures)
            {
                f.get();
            }

            return std::async(std::launch::deferred, [result]() { return result; });
        }

        std::future<Matrix<T>> multiply_async(const Matrix<T>& other, unsigned int blocks) const
        {
            if (M != other.N)
            {
                throw std::invalid_argument("Количество столбцов первой матрицы должно совпадать с количеством строк второй матрицы");
            }

            Matrix<T> result(N, other.M);
            auto multiply_blocks = [&](int start_row, int end_row) {
                for (int i = start_row; i < end_row; ++i)
                {
                    for (int j = 0; j < other.M; ++j)
                    {
                        T sum = 0;
                        for (int k = 0; k < M; ++k)
                        {
                            sum += this->get(i, k) * other.get(k, j);
                        }
                        result.set(i, j, sum);
                    }
                }
            };

            std::vector<std::future<void>> futures;
            for (int start_row = 0; start_row < N; start_row += static_cast<int>(blocks))
            {
                int end_row = std::min(start_row + static_cast<int>(blocks), N);
                futures.push_back(std::async(std::launch::async, multiply_blocks, start_row, end_row));
            }

            for (auto& f : futures)
            {
                f.get();
            }

            return std::async(std::launch::deferred, [result]() { return result; });
        }
};

template<typename T>
std::ostream &operator<<(std::ostream &outresult,const Matrix<T> &matr)
{
    const int N=matr.getN();
    const int M=matr.getM();

    outresult << '[' << std::endl;
    for (int indexN=0;indexN<N;indexN++)
    {
        outresult << " [ ";
        for (int indexM=0;indexM<M-1;indexM++)
        {
            outresult << matr.get(indexN,indexM) << ", ";
        }
        outresult << matr.get(indexN,M-1) << " ";
        outresult << "]" << std::endl;
    }
    outresult << ']' << std::endl;

    return outresult;
}

template <typename T>
std::istream &operator>>(std::istream &inmatrix,Matrix<T> &matr)
{
    int newN,newM;
    inmatrix >> newN >> newM;
    matr.changeSize(newN,newM);
    for (int indexN=0;indexN<newN;indexN++)
    {
        for (int indexM=0;indexM<newM;indexM++)
        {
            T elem;
            inmatrix >> elem;
            matr.set(indexN,indexM, elem);
        }
    }
    return inmatrix;
}

#endif  // MATRIX_H_
