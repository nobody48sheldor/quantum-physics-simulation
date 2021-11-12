#include <iostream>
#include <cmath>
#include <vector>
// #include <complex> //

using namespace std;

vector<vector<double>> psiinit(int L, int n, double alpha){
    vector<vector<double>> yj;
    for (int j = 0; j < n; j++)
    {
        vector<double> xi;
        for (int i = 0; i < n; i++)
        {
            double value = exp(-(pow((i-(L/4)), 2) + (pow((j-(L/4)), 2)))/alpha) / (sqrt(2)*3.14159*alpha);
            xi.push_back(value);
        }
        yj.push_back(xi);
    }
    return yj;
}

vector<vector<vector<double>>> psi_func(int L, double alpha, double m, double hbar, int n, int nt, vector<vector<double>> psi_init)
{
    double dx = 2*L/n;
    double dy = 2*L/n;
    double dt = 2*L/nt;

    int V = 0;

    double value;

    vector<vector<vector<double>>> psi;
    psi.push_back(psi_init);

    for (int k = 1; k < nt-2; k++)
    {
        vector<vector<double>> psi_y;
        psi_y.push_back(psi[k][0]);

        for (int j = 1; j < n-1; j++)
        {
            vector<double> psi_x;
            psi_x.push_back(psi[k][j][0]);

            for (int i = 1; i < n-1; i++)
            {
                value = psi[k-1][j][i] + (1/hbar) * ( ((hbar*hbar/(2*m))*(psi[k][j][i+1] - 2*psi[k][j][i] + psi[i][j][i-1])/(dx*dx) + (psi[k][j+1][i] - 2*psi[k][j][i] + psi[k][j-1][i])/(dy*dy)) + V*psi[k][j][i])*dt;
                psi_x.push_back(value);
            }
            psi_x.push_back(value);
            psi_y.push_back(psi_x);
        }
    psi.push_back(psi_y);
    }
    return psi;
}

int main(){
    int L = 10;
    int n = 400;
    int nt = 200*n;
    double alpha = 1;
    double m = 1;
    double hbar = 1;

    double x[n] = {};
    double y[n] = {};
    double t[nt] = {};

    vector<vector<double>> psi_init = psiinit(L, n, alpha);
    vector<vector<vector<double>>> psi = psi_func(L, alpha, m, hbar, n, nt, psi_init);

    cout << "done" <<endl;
    return 0;
}
