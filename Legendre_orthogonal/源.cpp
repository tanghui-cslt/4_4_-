#include <iostream>
#include <iomanip>
using namespace std;

/*****************��ҵҪ��****************

***  ����ָ������ex�����õ�����3�ζ���ʽ�ı����ʽ
***  ֱ�Ӹ���������x=1/2��ֵ��֤���

****************************************/

/*
input :
		(1) ���½�
		(2) 1/2
output:
		(1) ���õ���������ʽ�ı����ʽ
		(2) 1/2�Ľ��ӽ��
*/


void scanf_data(double upper_lower[], double x);
void solution(double upper_lower[], double x);
void calc_diag_mat(double upper_lower[], double diag_mat[]);

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
	cout << "\n\n********������������½�:**********\n";
	for (size_t i = 0; i < 2; i++)
	{
		cin >> upper_lower[i];
	}

	cout << "\n\n********��������֤����:**********\n";
	cin >> x;
}

