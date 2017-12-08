#include <iostream>
#include <iomanip>
using namespace std;

/*****************作业要求****************

***  构造指数函数ex的勒让德正交3次多项式的表达形式
***  直接给出，代入x=1/2的值验证结果

****************************************/

/*
input :
		(1) 上下界
		(2) 1/2
output:
		(1) 勒让德正交多项式的表达形式
		(2) 1/2的近视结果
*/


void scanf_data(double upper_lower[], double x);
void solution(double upper_lower[], double x);
void calc_diag_mat(double upper_lower[], double diag_mat[]);

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

	calc_diag_mat(upper_lower, diag_mat);



	for (int i = 0; i < 4; i++)
	{
		cout<<diag_mat[i]<<" ";
	}
}


void calc_diag_mat(double upper_lower[], double diag_mat[])
{
	for (int  i = 0; i < 4; i++)
	{
		diag_mat[i] = 2.0 / (2 * i + 1);
	}
}
void scanf_data(double upper_lower[], double x)
{
	cout << "\n\n********请输入积分上下界:**********\n";
	for (size_t i = 0; i < 2; i++)
	{
		cin >> upper_lower[i];
	}

	cout << "\n\n********请输入验证数据:**********\n";
	cin >> x;
}

