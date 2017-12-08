#include <iostream>
#include <iomanip>
using namespace std;

/*****************作业要求****************

***  构造ex的勒让德正交n次多项式的表达形式,用递归计算结果
***  直接给出，代入x=1/2的值验证结果

****************************************/

/*
input :
	(1) 上下界 (-1，1)
	(2) 迭代次数 10 （注意次数太高递归需要的时间较长 ）
	(3) 测试数据 0.5
output:
	(1) 1/2的近视结果
*/

void scanf_data(double upper_lower[], double &x, int &n);		//输入数据
void solution(double upper_lower[], double x, int n);			//求解
void calc_diag_mat(double upper_lower[], double diag_mat[], int n);//计算对角元素的值
void calc_a(double diag_mat[], double a[], double b[], int n);	//多项式系数
void calc_b(double upper_lower[], double b[], int n);			//计算右端的值
double calc_x(double a[], double x, int n);						//计算x的近视值
double integral_f(double temp_a[], int n);						//计算积分
double Pn_recusive(double x, int n);							//递归求解

int main()
{
	double x = 0;								//验证值
	double upper_lower[2] = { 0 };				//上下界
	int n = 0;									//多项式的阶数
	scanf_data(upper_lower, x, n);				//输入数据

	solution(upper_lower, x, n);				//求解

	getchar();
	getchar();

	return 0;
}

void solution(double upper_lower[], double x, int n)
{
	double *diag_mat = new double[n];			//左端对角阵
	double *b = new double[n];					//右端项
	double *a = new double[n];					//多项式系数

	calc_diag_mat(upper_lower, diag_mat, n);

	calc_b(upper_lower, b, n);					//使用函数指针计算右端项

	calc_a(diag_mat, a, b, n);					//计算多项式系数

	cout << "\n\n*******勒让德正交"<<n<<"次多项式的结果**********\n";
	cout << "\t" << calc_x(a, x, n);
}

void calc_b(double upper_lower[], double b[], int n)
{
	for (int i = 0; i < n; i++)
	{
		b[i] = integral_f(upper_lower, i);
	}

}

double calc_x(double a[], double x, int n)
{
	double res = 0.0;
	for (int i = 0; i < n; i++)
	{
		res += a[i] * Pn_recusive(x, i);
	}

	return res;
}

void calc_diag_mat(double upper_lower[], double diag_mat[], int n)
{
	for (int i = 0; i < n; i++)
	{
		diag_mat[i] = 2.0 / (2 * i + 1);
	}
}
void calc_a(double diag_mat[], double a[], double b[], int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] = b[i] / diag_mat[i];
	}
}

double Pn_recusive(double x, int n)
{
	double t = (double)n - 1;

	if (n == 0)
	{
		return 1;
	}
	else if (n == 1)
	{
		return x;
	}
	else
	{
		return (2.0*t + 1.0) / (t*1.0 + 1.0) *x *Pn_recusive(x, n - 1) - (t*1.0) / (t*1.0 + 1.0)*Pn_recusive(x, n - 2);
	}
}

double integral_f(double temp_a[], int n)
{
	double piecewise = 1000;			//分段个数
	double len = temp_a[1] - temp_a[0];	//区间的长度
	double sum = 0;						//积分的结果
	double piece_len = len / piecewise;	//每个小段的长度

	double x = temp_a[0];				//梯形的上底的x
	double x_1 = temp_a[0];				//梯形的下底的x
	double fx = 0;						//梯形的上底
	double(*fun)(double x, int n);		//函数指针
	fun = &Pn_recusive;

	double fx_1 = (*fun)(x_1, n)*exp(x_1);//梯形的下底
	for (int i = 1; i <= piecewise; i++)
	{
		x = x_1;
		x_1 = i * piece_len + temp_a[0];
		fx = fx_1;
		fx_1 = (*fun)(x_1, n)*exp(x_1);

		sum += 0.5*piece_len*(fx + fx_1);
	}
	return sum;
}


void scanf_data(double upper_lower[], double &x, int &n)
{
	cout << "\n\n********请输入积分上下界:**********\n";
	for (size_t i = 0; i < 2; i++)
	{
		cin >> upper_lower[i];
	}

	cout << "\n\n********请输入多项式阶数:**********\n";
	cin >> n;

	cout << "\n\n********请输入验证数据:**********\n";
	cin >> x;

}

