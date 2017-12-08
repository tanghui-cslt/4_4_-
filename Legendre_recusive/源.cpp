#include <iostream>
#include <iomanip>
using namespace std;

/*****************��ҵҪ��****************

***  ����ex�����õ�����n�ζ���ʽ�ı����ʽ,�õݹ������
***  ֱ�Ӹ���������x=1/2��ֵ��֤���

****************************************/

/*
input :
	(1) ���½� (-1��1)
	(2) �������� 10 ��ע�����̫�ߵݹ���Ҫ��ʱ��ϳ� ��
	(3) �������� 0.5
output:
	(1) 1/2�Ľ��ӽ��
*/

void scanf_data(double upper_lower[], double &x, int &n);		//��������
void solution(double upper_lower[], double x, int n);			//���
void calc_diag_mat(double upper_lower[], double diag_mat[], int n);//����Խ�Ԫ�ص�ֵ
void calc_a(double diag_mat[], double a[], double b[], int n);	//����ʽϵ��
void calc_b(double upper_lower[], double b[], int n);			//�����Ҷ˵�ֵ
double calc_x(double a[], double x, int n);						//����x�Ľ���ֵ
double integral_f(double temp_a[], int n);						//�������
double Pn_recusive(double x, int n);							//�ݹ����

int main()
{
	double x = 0;								//��ֵ֤
	double upper_lower[2] = { 0 };				//���½�
	int n = 0;									//����ʽ�Ľ���
	scanf_data(upper_lower, x, n);				//��������

	solution(upper_lower, x, n);				//���

	getchar();
	getchar();

	return 0;
}

void solution(double upper_lower[], double x, int n)
{
	double *diag_mat = new double[n];			//��˶Խ���
	double *b = new double[n];					//�Ҷ���
	double *a = new double[n];					//����ʽϵ��

	calc_diag_mat(upper_lower, diag_mat, n);

	calc_b(upper_lower, b, n);					//ʹ�ú���ָ������Ҷ���

	calc_a(diag_mat, a, b, n);					//�������ʽϵ��

	cout << "\n\n*******���õ�����"<<n<<"�ζ���ʽ�Ľ��**********\n";
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
	double piecewise = 1000;			//�ֶθ���
	double len = temp_a[1] - temp_a[0];	//����ĳ���
	double sum = 0;						//���ֵĽ��
	double piece_len = len / piecewise;	//ÿ��С�εĳ���

	double x = temp_a[0];				//���ε��ϵ׵�x
	double x_1 = temp_a[0];				//���ε��µ׵�x
	double fx = 0;						//���ε��ϵ�
	double(*fun)(double x, int n);		//����ָ��
	fun = &Pn_recusive;

	double fx_1 = (*fun)(x_1, n)*exp(x_1);//���ε��µ�
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
	cout << "\n\n********������������½�:**********\n";
	for (size_t i = 0; i < 2; i++)
	{
		cin >> upper_lower[i];
	}

	cout << "\n\n********���������ʽ����:**********\n";
	cin >> n;

	cout << "\n\n********��������֤����:**********\n";
	cin >> x;

}

