#include <iostream>
#include <iomanip>
using namespace std;
#define Error 1e-12

void scanf_data(double temp_a[]);
void scanf_data(double temp_a[], int &n);
void successive_half(double temp_a[], int n);
double f_x(double x);
double* integral_f(double temp_a[], int n, double(*fun)(double));

int main()
{
	double temp_a[2] = { 0 };				// ���½�
	int n = 2;								// ��������
	scanf_data(temp_a);
	successive_half(temp_a, n);
	getchar();
	getchar();
}


void successive_half(double temp_a[], int n)
{
	double error = 100;				//��ʼ���ֵ ��Ҫ >0
	double T_n = 0;
	double T_2n = 0;
	double S_2n = 0;
	double H_n = 0;
	double(*fun)(double) = f_x;

	double *temp = NULL;

	temp = integral_f(temp_a, n, fun);
	T_2n = temp[0];
	cout << "******�ֱ���n, Tn,T2n,S2n******\n";
	while (fabs(error) > Error)
	{
		T_n = T_2n;
		temp = integral_f(temp_a, n, fun);

		T_2n = 0.5*temp[0] + 0.5*temp[1];

		S_2n = 4.0 / 3.0 * T_2n - T_n / 3.0;

		error = T_2n - T_n;
		//T_n = T_2n;
		n = n * 2;
		cout << "n = " << n << " ";
		cout <<setprecision(14) << " "<<T_n <<" "<<T_2n<<" "<<S_2n<<endl;
		//getchar();

	}

	cout << "\n\n************���Ϊ*********\n";
	cout  << "T2n = " << T_2n << endl;
	cout  << "S2n = " << S_2n << endl;
}

double f_x(double x)
{
	return exp(x)*cos(x);
}

double* integral_f(double temp_a[], int n, double(*fun)(double))
{
	double piecewise = n;				//�ֶθ���
	double len = temp_a[1] - temp_a[0];	//����ĳ���
	double sum = 0;						//���ֵĽ��
	double piece_len = len / piecewise;	//ÿ��С�εĳ���

	double x = temp_a[0];				//���ε��ϵ׵�x
	double x_1 = temp_a[0];				//���ε��µ׵�x
	double fx = 0;						//���ε��ϵ�
	double fx_1 = 1.0 / (1 + x_1*x_1);	//���ε��µ�

	double fx_0_1 = 0;					//�е�ĺ���ֵ
	double rec_area = 0;				//���ε����

	double *area = new double[2];
	for (int i = 1; i <= piecewise; i++)
	{
		x = x_1;
		x_1 = i * piece_len + temp_a[0];
		fx = fx_1;
		fx_1 = (*fun)(x_1);

		sum += 0.5*piece_len*(fx + fx_1);
		//cout << "x = " << 0.5*(x + x_1) << endl;
		//getchar();
		fx_0_1 = (*fun)(0.5*(x + x_1));
		rec_area += fx_0_1 * piece_len;
	}
	//cout << "\n\n************���Ϊ*********\n";
	//cout << sum <<" rec = "<< rec_area << endl;

	area[0] = sum;
	area[1] = rec_area;
	return area;
}

void  scanf_data(double temp_a[])
{
	cout << "\n\n********������������½�:**********\n";

	for (size_t i = 0; i < 2; i++)
	{
		cin >> temp_a[i];
	}
}
void scanf_data(double temp_a[], int &n)
{
	cout << "\n\n********������������½�:**********\n";

	for (size_t i = 0; i < 2; i++)
	{
		cin >> temp_a[i];
	}

	cout << "\n\n********�������ʼ��������:**********\n";

	cin >> n;
}
