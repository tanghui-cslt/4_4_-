#include <iostream>
#include <iomanip>
using namespace std;

/*****************作业要求****************

***  构造指数函数ex的勒让德正交3次多项式的表达形式
***  直接给出，代入x=1/2的值验证结果

****************************************/

/*
input :
		(1) 上下界 (-1，1)
		(2) 测试数据 0.5
output:
		(1) 1/2的近视结果
*/

void scanf_data(double upper_lower[], double &x);				//输入数据
void solution(double upper_lower[], double x);					//求解
void calc_diag_mat(double upper_lower[], double diag_mat[]);	//计算对角元素的值
void calc_a(double diag_mat[], double a[], double b[], int n);	//多项式系数
void calc_b(double upper_lower[], double b[]);					//计算右端的值
double calc_x(double a[], double x, int n);
double integral_f(double temp_a[], double(*fun)(double));		//计算积分
double P0(double x);
double P1(double x);
double P2(double x);
double P3(double x);

int main()
{
	double x = 0;								//验证值
	double upper_lower[2] = { 0 };				//上下界

	scanf_data(upper_lower, x);					//输入数据

	solution(upper_lower, x);					//求解

	getchar();
	getchar();

	return 0;
}

void solution(double upper_lower[], double x)
{
	double diag_mat[4] = { 0 };					//左端对角阵
	double b[4] = { 0 };						//右端项
	double a[4] = { 0 };						//多项式系数
	int n = 4;									//计算维数

	calc_diag_mat(upper_lower, diag_mat);

	calc_b(upper_lower, b);						//使用函数指针计算右端项

	calc_a(diag_mat, a, b, n);					//计算多项式系数

	cout << "\n\n*******勒让德正交3次多项式的结果**********\n";
	cout << "\t" << calc_x(a, x, n);
}

void calc_b(double upper_lower[], double b[])
{
	double(*fun)(double);						//定义函数指针，计算右端项

	fun = &P0;
	b[0] = integral_f(upper_lower, fun);

	fun = &P1;
	b[1] = integral_f(upper_lower, fun);

	fun = &P2;
	b[2] = integral_f(upper_lower, fun);

	fun = &P3;
	b[3] = integral_f(upper_lower, fun);
}

double calc_x(double a[], double x, int n)
{
	double res = 0.0;

	res += a[0] * P0(x);
	res += a[1] * P1(x);
	res += a[2] * P2(x);
	res += a[3] * P3(x);

	return res;
}

void calc_diag_mat(double upper_lower[], double diag_mat[])
{
	for (int i = 0; i < 4; i++)
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
double P0(double x)
{
	return 1;
}
double P1(double x)
{
	return x;
}
double P2(double x)
{
	return (3 * x*x - 1.0) / 2.0;
}
double P3(double x)
{
	return (5 * pow(x, 3) - 3 * x) / 2.0;
}

double integral_f(double temp_a[], double(*fun)(double))
{
	double piecewise = 4000;			//分段个数
	double len = temp_a[1] - temp_a[0];	//区间的长度
	double sum = 0;						//积分的结果
	double piece_len = len / piecewise;	//每个小段的长度

	double x = temp_a[0];				//梯形的上底的x
	double x_1 = temp_a[0];				//梯形的下底的x
	double fx = 0;						//梯形的上底
	double fx_1 = (*fun)(x_1)*exp(x_1);	//梯形的下底
	for (int i = 1; i <= piecewise; i++)
	{
		x = x_1;
		x_1 = i * piece_len + temp_a[0];
		fx = fx_1;
		fx_1 = (*fun)(x_1)*exp(x_1);

		sum += 0.5*piece_len*(fx + fx_1);
	}
	//cout << "\n\n************结果为*********\n";
	//cout << sum << endl;
	return sum;
}


void scanf_data(double upper_lower[], double &x)
{
	cout << "\n\n********请输入积分上下界:**********\n";
	for (size_t i = 0; i < 2; i++)
	{
		cin >> upper_lower[i];
	}

	cout << "\n\n********请输入验证数据:**********\n";
	cin >> x;
}

