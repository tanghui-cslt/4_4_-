#include <iostream>
#include <iomanip>
using namespace std;

/*****************��ҵҪ��****************

***  ����ָ������ex�����õ�����3�ζ���ʽ�ı����ʽ
***  ֱ�Ӹ���������x=1/2��ֵ��֤���

****************************************/

/*
input :
		(1) ���½� (-1��1)
		(2) �������� 0.5
output:
		(1) 1/2�Ľ��ӽ��
*/

void scanf_data(double upper_lower[], double &x);				//��������
void solution(double upper_lower[], double x);					//���
void calc_diag_mat(double upper_lower[], double diag_mat[]);	//����Խ�Ԫ�ص�ֵ
void calc_a(double diag_mat[], double a[], double b[], int n);	//����ʽϵ��
void calc_b(double upper_lower[], double b[]);					//�����Ҷ˵�ֵ
double calc_x(double a[], double x, int n);
double integral_f(double temp_a[], double(*fun)(double));		//�������
double P0(double x);
double P1(double x);
double P2(double x);
double P3(double x);

int main()
{
	double x = 0;								//��ֵ֤
	double upper_lower[2] = { 0 };				//���½�

	scanf_data(upper_lower, x);					//��������

	solution(upper_lower, x);					//���

	getchar();
	getchar();

	return 0;
}

void solution(double upper_lower[], double x)
{
	double diag_mat[4] = { 0 };					//��˶Խ���
	double b[4] = { 0 };						//�Ҷ���
	double a[4] = { 0 };						//����ʽϵ��
	int n = 4;									//����ά��

	calc_diag_mat(upper_lower, diag_mat);

	calc_b(upper_lower, b);						//ʹ�ú���ָ������Ҷ���

	calc_a(diag_mat, a, b, n);					//�������ʽϵ��

	cout << "\n\n*******���õ�����3�ζ���ʽ�Ľ��**********\n";
	cout << "\t" << calc_x(a, x, n);
}

void calc_b(double upper_lower[], double b[])
{
	double(*fun)(double);						//���庯��ָ�룬�����Ҷ���

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
	double piecewise = 4000;			//�ֶθ���
	double len = temp_a[1] - temp_a[0];	//����ĳ���
	double sum = 0;						//���ֵĽ��
	double piece_len = len / piecewise;	//ÿ��С�εĳ���

	double x = temp_a[0];				//���ε��ϵ׵�x
	double x_1 = temp_a[0];				//���ε��µ׵�x
	double fx = 0;						//���ε��ϵ�
	double fx_1 = (*fun)(x_1)*exp(x_1);	//���ε��µ�
	for (int i = 1; i <= piecewise; i++)
	{
		x = x_1;
		x_1 = i * piece_len + temp_a[0];
		fx = fx_1;
		fx_1 = (*fun)(x_1)*exp(x_1);

		sum += 0.5*piece_len*(fx + fx_1);
	}
	//cout << "\n\n************���Ϊ*********\n";
	//cout << sum << endl;
	return sum;
}


void scanf_data(double upper_lower[], double &x)
{
	cout << "\n\n********������������½�:**********\n";
	for (size_t i = 0; i < 2; i++)
	{
		cin >> upper_lower[i];
	}

	cout << "\n\n********��������֤����:**********\n";
	cin >> x;
}

