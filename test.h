#include <iostream>
#include <vector>
#include <cmath> // 为了使用 abs() 函数

// 打印矩阵的辅助函数
void print_matrix(const std::vector<std::vector<double>>& mat) {
    for (const auto& row : mat) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

// 高斯消元法解线性方程组 A * x = b
std::vector<double> gaussian_elimination(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int n = A.size();

    // 将 A 和 b 合并成一个增广矩阵
    for (int i = 0; i < n; ++i) {
        A[i].push_back(b[i]);
    }

    // 执行高斯消元
    for (int i = 0; i < n; ++i) {
        // 寻找当前列最大值作为主元素
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                max_row = k;
            }
        }

        // 交换行
        std::swap(A[i], A[max_row]);

        // 归一化当前行
        double pivot = A[i][i];
        for (int j = i; j < n + 1; ++j) {
            A[i][j] /= pivot;
        }

        // 消去其他行
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = i; j < n + 1; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
            }
        }
    }

    // 从增广矩阵中提取解
    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = A[i][n]; // 最后一列是解向量
    }
    return result;
}

// 封装的多项式拟合函数
std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y, int degree) {
    int n = x.size();

    // 构造矩阵 A 和向量 b
    std::vector<std::vector<double>> A(n, std::vector<double>(degree + 1)); // A 是 n x (degree+1)
    std::vector<double> b(n);

    // 构建矩阵 A 和向量 b
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= degree; ++j) {
            A[i][j] = std::pow(x[i], j);  // x[i] 的不同幂次
        }
        b[i] = y[i];
    }

    // 使用高斯消元法解方程 A * coeff = b
    return gaussian_elimination(A, b);
}

int main() {
    // 输入数据点
    std::vector<double> x = {0, 3, 6};
    std::vector<double> y = {400, 400, 400};

    // 调用 polyfit 函数拟合 2 阶多项式
    std::vector<double> coeff = polyfit(x, y, 2);

    // 输出系数 a, b, c
    std::cout << "a: " << coeff[2] << ", b: " << coeff[1] << ", c: " << coeff[0] << std::endl;

    return 0;
}


