#include <iostream>
#include <cmath>
#include <vector>
// #include <complex> //

using namespace std;

auto init(int n, double l, double alpha)
{
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double value = exp(-(pow((i-(l/4)), 2) + (pow((j-(l/4)), 2)))/alpha) / (sqrt(2)*3.14159*alpha);
            xarray[i] = value;
        }
        yarray[j] = xarray;
    }
    return yarray;
}

auto psi_func(double l, double alpha, double m, double hbar, int n, int nt, vector<vector<double>> psi_init, int treshold)
{
    double dx = 2*l/n;
    double dy = 2*l/n;
    double dt = 2*l/nt;

    int v = 0;

    double value;

    vector<vector<vector<double>>> psi;
    psi.resize(nt);
    vector<vector<double>> psi_y;
    psi_y.resize(n);
    vector<double> psi_x;
    psi_x.resize(n);


    vector<vector<double>> empty;
    empty.resize(1);
    vector<double> empty_;
    empty_.resize(1);

    empty_[0] = 0;
    empty[0] = empty_;

    psi[0] = psi_init;

    for (int k = 1; k < nt; k++)
    {
        if (k > treshold)
        {
            psi[k - treshold + 1] = empty;
        }
        cout << k << " / " << nt << endl;

        psi_y[0] = psi[k-1][0];

        for (int j = 1; j < n-1; j++)
        {
            psi_x[0] = psi[k-1][j][0];

            for (int i = 1; i < n-1; i++)
            {
                value = psi[k-1][j][i] + ((1/hbar) * ( ((hbar*hbar/(2*m))*(psi[k-1][j][i+1] - 2*psi[k-1][j][i] + psi[k-1][j][i-1])/(dx*dx) + (psi[k-1][j+1][i] - 2*psi[k-1][j][i] + psi[k-1][j-1][i])/(dy*dy)) + v*psi[k-1][j][i])*dt);
                psi_x[i] = value;
            }
            psi_x[n-1] = value;
            psi_y[j] = psi_x;
        }
        psi_y[n-1] = psi_x;
        psi[k] = psi_y;
    }
    return psi;
}

int main()
{
    double l = 10;
    const int n = 800;
    const int nt = 15000;
    double alpha = 1;
    double m = 1;
    double hbar = 1;

    int treshold = 20;

    vector<vector<double>> psi_init = init(n, l, alpha);
    vector<vector<vector<double>>> psi = psi_func(l, alpha, m, hbar, n, nt, psi_init, treshold);
    psi.resize(nt);
    cout << "done" <<endl;

    return 0;
}
