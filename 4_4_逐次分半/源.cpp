#include <iostream>
#include <iomanip>
using namespace std;
#define Error 1e-12

void scanf_data(double temp_a[]);			//重载函数
void scanf_data(double temp_a[], int &n);	//输入数据
void successive_half(double temp_a[], int n);
double f_x(double x);
double* integral_f(double temp_a[], int n, double(*fun)(double));

int main()
{
	double temp_a[2] = { 0 };				// 上下界
	int n = 2;								// 迭代次数
	scanf_data(temp_a);
	successive_half(temp_a, n);
	getchar();
	getchar();
}


void successive_half(double temp_a[], int n)
{
	double error = 100;				//初始误差值 需要 >0
	double T_n = 0;
	double T_2n = 0;
	double S_2n = 0;
	double H_n = 0;
	double S_n = 0,temp_s_n;
	double(*fun)(double) = f_x;
	double times = 0, S_times = 0;
	double *temp = NULL;

	temp = integral_f(temp_a, n, fun);
	T_2n = temp[0];
	cout << "\n\n\ttimes\tn\t\tT2n\t\tS2n\n";
	while (fabs(T_2n - T_n) > Error)
	{
		times++;
		T_n = T_2n;
		S_n = S_2n;
		temp = integral_f(temp_a, n, fun);

		T_2n = 0.5*temp[0] + 0.5*temp[1];
		S_2n = 4.0 / 3.0 * T_2n - T_n / 3.0;

		cout << "\t" << times << "\t" << n;
		cout << setprecision(14) << "\t" << T_2n << "\t";
		if (fabs(S_2n - S_n) > Error)
		{
			S_times = times;
			temp_s_n = S_2n;
			cout << S_2n << endl;
		}
		else	cout << endl;

		n = n * 2;
	}

	cout << "\n\n************结果为*********\n";
	cout << "\tT_2n运行" << times << "次，结果为 = " << T_2n << endl;
	cout << "\tS_2n运行" << S_times << "次，结果为 = " <<temp_s_n << endl;
}

double f_x(double x)
{
	return exp(x)*cos(x);
}

double* integral_f(double temp_a[], int n, double(*fun)(double))
{
	double piecewise = n;				//分段个数
	double len = temp_a[1] - temp_a[0];	//区间的长度
	double sum = 0;						//积分的结果
	double piece_len = len / piecewise;	//每个小段的长度

	double x = temp_a[0];				//梯形的上底的x
	double x_1 = temp_a[0];				//梯形的下底的x
	double fx = 0;						//梯形的上底
	double fx_1 = (*fun)(x_1);			//梯形的下底

	double fx_0_1 = 0;					//中点的函数值
	double rec_area = 0;				//矩形的面积

	double *area = new double[2];
	for (int i = 1; i <= piecewise; i++)
	{
		x = x_1;
		x_1 = i * piece_len + temp_a[0];
		fx = fx_1;
		fx_1 = (*fun)(x_1);

		sum += 0.5*piece_len*(fx + fx_1);
		fx_0_1 = (*fun)(0.5*(x + x_1));
		rec_area += fx_0_1 * piece_len;
	}

	area[0] = sum;
	area[1] = rec_area;
	return area;
}

void  scanf_data(double temp_a[])
{
	cout << "\n\n********请输入积分上下界:**********\n";

	for (size_t i = 0; i < 2; i++)
	{
		cin >> temp_a[i];
	}
}
void scanf_data(double temp_a[], int &n)
{
	cout << "\n\n********请输入积分上下界:**********\n";

	for (size_t i = 0; i < 2; i++)
	{
		cin >> temp_a[i];
	}

	cout << "\n\n********请输入初始迭代次数:**********\n";

	cin >> n;
}
