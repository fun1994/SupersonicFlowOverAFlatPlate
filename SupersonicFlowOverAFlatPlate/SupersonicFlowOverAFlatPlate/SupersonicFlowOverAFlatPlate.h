#pragma once
#include <iostream>
#include "Data.h"

class SupersonicFlowOverAFlatPlate {
	double gamma;
	double R;
	double c_p;
	double c_v;
	double mu0;
	double T0;
	double Pr;
	double Ma_inf;
	double p_inf;
	double T_inf;
	double T_w;
	std::string wall;
	int Nx;
	int Ny;
	double dx;
	double dy;
	double K;
	double tol;
	void initialize(Data& data, std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& e, std::vector<std::vector<double>>& mu, std::vector<std::vector<double>>& k, std::vector<std::vector<double>>& Ma, std::vector<std::vector<double>>& U1, std::vector<std::vector<double>>& U2, std::vector<std::vector<double>>& U3, std::vector<std::vector<double>>& U5, std::vector<std::vector<double>>& tau_xx0, std::vector<std::vector<double>>& tau_xx1, std::vector<std::vector<double>>& tau_xy0, std::vector<std::vector<double>>& tau_xy1, std::vector<std::vector<double>>& tau_xy2, std::vector<std::vector<double>>& tau_xy3, std::vector<std::vector<double>>& tau_yy0, std::vector<std::vector<double>>& tau_yy1, std::vector<std::vector<double>>& q_x0, std::vector<std::vector<double>>& q_x1, std::vector<std::vector<double>>& q_y0, std::vector<std::vector<double>>& q_y1, std::vector<std::vector<double>>& E1, std::vector<std::vector<double>>& E2, std::vector<std::vector<double>>& E3, std::vector<std::vector<double>>& E5, std::vector<std::vector<double>>& F1, std::vector<std::vector<double>>& F2, std::vector<std::vector<double>>& F3, std::vector<std::vector<double>>& F5, std::vector<std::vector<double>>& dU1dt, std::vector<std::vector<double>>& dU2dt, std::vector<std::vector<double>>& dU3dt, std::vector<std::vector<double>>& dU5dt, std::vector<std::vector<double>>& U1_bar, std::vector<std::vector<double>>& U2_bar, std::vector<std::vector<double>>& U3_bar, std::vector<std::vector<double>>& U5_bar, std::vector<std::vector<double>>& rho_bar, std::vector<std::vector<double>>& u_bar, std::vector<std::vector<double>>& v_bar, std::vector<std::vector<double>>& T_bar, std::vector<std::vector<double>>& p_bar, std::vector<std::vector<double>>& e_bar, std::vector<std::vector<double>>& mu_bar, std::vector<std::vector<double>>& k_bar, std::vector<std::vector<double>>& Ma_bar, std::vector<std::vector<double>>& E1_bar, std::vector<std::vector<double>>& E2_bar, std::vector<std::vector<double>>& E3_bar, std::vector<std::vector<double>>& E5_bar, std::vector<std::vector<double>>& F1_bar, std::vector<std::vector<double>>& F2_bar, std::vector<std::vector<double>>& F3_bar, std::vector<std::vector<double>>& F5_bar, std::vector<std::vector<double>>& dU1dt_bar, std::vector<std::vector<double>>& dU2dt_bar, std::vector<std::vector<double>>& dU3dt_bar, std::vector<std::vector<double>>& dU5dt_bar, std::vector<std::vector<double>>& dU1dt_av, std::vector<std::vector<double>>& dU2dt_av, std::vector<std::vector<double>>& dU3dt_av, std::vector<std::vector<double>>& dU5dt_av) {
		data.x.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			data.x[i] = i * dx;
		}
		data.y.resize(Ny + 1);
		for (int i = 0; i < Ny + 1; i++) {
			data.y[i] = i * dy;
		}
		u.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			u[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				u[i][j] = Ma_inf * sqrt(gamma * R * T_inf);
			}
		}
		v.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			v[i].resize(Ny + 1);
		}
		T.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			T[i].resize(Ny + 1);
		}
		T[0][0] = T_inf;
		for (int i = 1; i < Nx + 1; i++) {
			T[i][0] = T_w;
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				T[i][j] = T_inf;
			}
		}
		p.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			p[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				p[i][j] = p_inf;
			}
		}
		rho.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			rho[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				rho[i][j] = p[i][j] / (R * T[i][j]);
			}
		}
		e.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			e[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				e[i][j] = c_v * T[i][j];
			}
		}
		mu.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			mu[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				mu[i][j] = mu0 * pow(T[i][j] / T0, 1.5) * (T0 + 110) / (T[i][j] + 110);
			}
		}
		k.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			k[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				k[i][j] = mu[i][j] * c_p / Pr;
			}
		}
		Ma.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			Ma[i].resize(Ny + 1);
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				Ma[i][j] = Ma_inf;
			}
		}
		U1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U1[i].resize(Ny + 1);
		}
		U2.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U2[i].resize(Ny + 1);
		}
		U3.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U3[i].resize(Ny + 1);
		}
		U5.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U5[i].resize(Ny + 1);
		}
		tau_xx0.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_xx0[i].resize(Ny + 1);
		}
		tau_xx1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_xx1[i].resize(Ny + 1);
		}
		tau_xy0.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_xy0[i].resize(Ny + 1);
		}
		tau_xy1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_xy1[i].resize(Ny + 1);
		}
		tau_xy2.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_xy2[i].resize(Ny + 1);
		}
		tau_xy3.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_xy3[i].resize(Ny + 1);
		}
		tau_yy0.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_yy0[i].resize(Ny + 1);
		}
		tau_yy1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			tau_yy1[i].resize(Ny + 1);
		}
		q_x0.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			q_x0[i].resize(Ny + 1);
		}
		q_x1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			q_x1[i].resize(Ny + 1);
		}
		q_y0.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			q_y0[i].resize(Ny + 1);
		}
		q_y1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			q_y1[i].resize(Ny + 1);
		}
		E1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E1[i].resize(Ny + 1);
		}
		E2.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E2[i].resize(Ny + 1);
		}
		E3.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E3[i].resize(Ny + 1);
		}
		E5.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E5[i].resize(Ny + 1);
		}
		F1.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F1[i].resize(Ny + 1);
		}
		F2.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F2[i].resize(Ny + 1);
		}
		F3.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F3[i].resize(Ny + 1);
		}
		F5.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F5[i].resize(Ny + 1);
		}
		dU1dt.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU1dt[i].resize(Ny + 1);
		}
		dU2dt.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU2dt[i].resize(Ny + 1);
		}
		dU3dt.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU3dt[i].resize(Ny + 1);
		}
		dU5dt.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU5dt[i].resize(Ny + 1);
		}
		U1_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U1_bar[i].resize(Ny + 1);
		}
		U2_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U2_bar[i].resize(Ny + 1);
		}
		U3_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U3_bar[i].resize(Ny + 1);
		}
		U5_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			U5_bar[i].resize(Ny + 1);
		}
		rho_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			rho_bar[i].resize(Ny + 1);
		}
		u_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			u_bar[i].resize(Ny + 1);
		}
		v_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			v_bar[i].resize(Ny + 1);
		}
		T_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			T_bar[i].resize(Ny + 1);
		}
		p_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			p_bar[i].resize(Ny + 1);
		}
		e_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			e_bar[i].resize(Ny + 1);
		}
		mu_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			mu_bar[i].resize(Ny + 1);
		}
		k_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			k_bar[i].resize(Ny + 1);
		}
		Ma_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			Ma_bar[i].resize(Ny + 1);
		}
		E1_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E1_bar[i].resize(Ny + 1);
		}
		E2_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E2_bar[i].resize(Ny + 1);
		}
		E3_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E3_bar[i].resize(Ny + 1);
		}
		E5_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			E5_bar[i].resize(Ny + 1);
		}
		F1_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F1_bar[i].resize(Ny + 1);
		}
		F2_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F2_bar[i].resize(Ny + 1);
		}
		F3_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F3_bar[i].resize(Ny + 1);
		}
		F5_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			F5_bar[i].resize(Ny + 1);
		}
		dU1dt_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU1dt_bar[i].resize(Ny + 1);
		}
		dU2dt_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU2dt_bar[i].resize(Ny + 1);
		}
		dU3dt_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU3dt_bar[i].resize(Ny + 1);
		}
		dU5dt_bar.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU5dt_bar[i].resize(Ny + 1);
		}
		dU1dt_av.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU1dt_av[i].resize(Ny + 1);
		}
		dU2dt_av.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU2dt_av[i].resize(Ny + 1);
		}
		dU3dt_av.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU3dt_av[i].resize(Ny + 1);
		}
		dU5dt_av.resize(Nx + 1);
		for (int i = 0; i < Nx + 1; i++) {
			dU5dt_av[i].resize(Ny + 1);
		}
	}
	void save(Data& data, double dt, std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& Ma) {
		if (data.t.empty()) {
			data.t.push_back(0);
		}
		else {
			data.t.push_back(data.t[data.t.size() - 1] + dt);
		}
		data.rho.push_back(rho);
		data.u.push_back(u);
		data.v.push_back(v);
		data.T.push_back(T);
		data.p.push_back(p);
		data.Ma.push_back(Ma);
	}
	double timeStep(std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& mu) {
		double a = sqrt(gamma * R * T[1][1]);
		double nu_prime = 4.0 / 3.0 * (gamma * mu[1][1] / Pr) / rho[1][1];
		double dt_CFL = 1 / (abs(u[1][1]) / dx + abs(v[1][1]) / dy + a * sqrt(1 / pow(dx, 2) + 1 / pow(dy, 2)) + 2 * nu_prime * (1 / pow(dx, 2) + 1 / pow(dy, 2)));
		double dt = K * dt_CFL;
		for (int i = 1; i < Nx; i++) {
			for (int j = 1; j < Ny; j++) {
				a = sqrt(gamma * R * T[i][j]);
				nu_prime = 4.0 / 3.0 * (gamma * mu[i][j] / Pr) / rho[i][j];
				dt_CFL = 1 / (abs(u[i][j]) / dx + abs(v[i][j]) / dy + a * sqrt(1 / pow(dx, 2) + 1 / pow(dy, 2)) + 2 * nu_prime * (1 / pow(dx, 2) + 1 / pow(dy, 2)));
				dt = dt < K * dt_CFL ? dt : K * dt_CFL;
			}
		}
		return dt;
	}
	void MacCormack(double dt, std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& e, std::vector<std::vector<double>>& mu, std::vector<std::vector<double>>& k, std::vector<std::vector<double>>& Ma, std::vector<std::vector<double>>& U1, std::vector<std::vector<double>>& U2, std::vector<std::vector<double>>& U3, std::vector<std::vector<double>>& U5, std::vector<std::vector<double>>& tau_xx0, std::vector<std::vector<double>>& tau_xx1, std::vector<std::vector<double>>& tau_xy0, std::vector<std::vector<double>>& tau_xy1, std::vector<std::vector<double>>& tau_xy2, std::vector<std::vector<double>>& tau_xy3, std::vector<std::vector<double>>& tau_yy0, std::vector<std::vector<double>>& tau_yy1, std::vector<std::vector<double>>& q_x0, std::vector<std::vector<double>>& q_x1, std::vector<std::vector<double>>& q_y0, std::vector<std::vector<double>>& q_y1, std::vector<std::vector<double>>& E1, std::vector<std::vector<double>>& E2, std::vector<std::vector<double>>& E3, std::vector<std::vector<double>>& E5, std::vector<std::vector<double>>& F1, std::vector<std::vector<double>>& F2, std::vector<std::vector<double>>& F3, std::vector<std::vector<double>>& F5, std::vector<std::vector<double>>& dU1dt, std::vector<std::vector<double>>& dU2dt, std::vector<std::vector<double>>& dU3dt, std::vector<std::vector<double>>& dU5dt, std::vector<std::vector<double>>& U1_bar, std::vector<std::vector<double>>& U2_bar, std::vector<std::vector<double>>& U3_bar, std::vector<std::vector<double>>& U5_bar, std::vector<std::vector<double>>& rho_bar, std::vector<std::vector<double>>& u_bar, std::vector<std::vector<double>>& v_bar, std::vector<std::vector<double>>& T_bar, std::vector<std::vector<double>>& p_bar, std::vector<std::vector<double>>& e_bar, std::vector<std::vector<double>>& mu_bar, std::vector<std::vector<double>>& k_bar, std::vector<std::vector<double>>& Ma_bar, std::vector<std::vector<double>>& E1_bar, std::vector<std::vector<double>>& E2_bar, std::vector<std::vector<double>>& E3_bar, std::vector<std::vector<double>>& E5_bar, std::vector<std::vector<double>>& F1_bar, std::vector<std::vector<double>>& F2_bar, std::vector<std::vector<double>>& F3_bar, std::vector<std::vector<double>>& F5_bar, std::vector<std::vector<double>>& dU1dt_bar, std::vector<std::vector<double>>& dU2dt_bar, std::vector<std::vector<double>>& dU3dt_bar, std::vector<std::vector<double>>& dU5dt_bar, std::vector<std::vector<double>>& dU1dt_av, std::vector<std::vector<double>>& dU2dt_av, std::vector<std::vector<double>>& dU3dt_av, std::vector<std::vector<double>>& dU5dt_av) {
		getTauxx(tau_xx0, u, v, mu, 0);
		getTauxy(tau_xy0, u, v, mu, 0);
		getQx(q_x0, T, k, 0);
		getE(E1, E2, E3, E5, rho, u, v, p, e, tau_xx0, tau_xy0, q_x0);
		getTauxy(tau_xy1, u, v, mu, 1);
		getTauyy(tau_yy0, u, v, mu, 0);
		getQy(q_y0, T, k, 0);
		getF(F1, F2, F3, F5, rho, u, v, p, e, tau_xy1, tau_yy0, q_y0);
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				dU1dt[i][j] = -(E1[i + 1][j] - E1[i][j]) / dx - (F1[i][j + 1] - F1[i][j]) / dy;
			}
		}
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				dU2dt[i][j] = -(E2[i + 1][j] - E2[i][j]) / dx - (F2[i][j + 1] - F2[i][j]) / dy;
			}
		}
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				dU3dt[i][j] = -(E3[i + 1][j] - E3[i][j]) / dx - (F3[i][j + 1] - F3[i][j]) / dy;
			}
		}
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				dU5dt[i][j] = -(E5[i + 1][j] - E5[i][j]) / dx - (F5[i][j + 1] - F5[i][j]) / dy;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U1_bar[i][j] = U1[i][j] + dU1dt[i][j] * dt;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U2_bar[i][j] = U2[i][j] + dU2dt[i][j] * dt;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U3_bar[i][j] = U3[i][j] + dU3dt[i][j] * dt;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U5_bar[i][j] = U5[i][j] + dU5dt[i][j] * dt;
			}
		}
		getPrimitiveVariables(rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, mu_bar, k_bar, Ma_bar, U1_bar, U2_bar, U3_bar, U5_bar);
		boundaryCondition(rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, mu_bar, k_bar, Ma_bar, U1_bar, U2_bar, U3_bar, U5_bar);
		getTauxx(tau_xx1, u_bar, v_bar, mu_bar, 1);
		getTauxy(tau_xy2, u_bar, v_bar, mu_bar, 2);
		getQx(q_x1, T_bar, k_bar, 1);
		getE(E1_bar, E2_bar, E3_bar, E5_bar, rho_bar, u_bar, v_bar, p_bar, e_bar, tau_xx1, tau_xy2, q_x1);
		getTauxy(tau_xy3, u_bar, v_bar, mu_bar, 3);
		getTauyy(tau_yy1, u_bar, v_bar, mu_bar, 1);
		getQy(q_y1, T_bar, k_bar, 1);
		getF(F1_bar, F2_bar, F3_bar, F5_bar, rho_bar, u_bar, v_bar, p_bar, e_bar, tau_xy3, tau_yy1, q_y1);
		for (int i = 1; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				dU1dt_bar[i][j] = -(E1_bar[i][j] - E1_bar[i - 1][j]) / dx - (F1_bar[i][j] - F1_bar[i][j - 1]) / dy;
			}
		}
		for (int i = 1; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				dU2dt_bar[i][j] = -(E2_bar[i][j] - E2_bar[i - 1][j]) / dx - (F2_bar[i][j] - F2_bar[i][j - 1]) / dy;
			}
		}
		for (int i = 1; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				dU3dt_bar[i][j] = -(E3_bar[i][j] - E3_bar[i - 1][j]) / dx - (F3_bar[i][j] - F3_bar[i][j - 1]) / dy;
			}
		}
		for (int i = 1; i < Nx + 1; i++) {
			for (int j = 1; j < Ny + 1; j++) {
				dU5dt_bar[i][j] = -(E5_bar[i][j] - E5_bar[i - 1][j]) / dx - (F5_bar[i][j] - F5_bar[i][j - 1]) / dy;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				dU1dt_av[i][j] = (dU1dt[i][j] + dU1dt_bar[i][j]) / 2;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				dU2dt_av[i][j] = (dU2dt[i][j] + dU2dt_bar[i][j]) / 2;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				dU3dt_av[i][j] = (dU3dt[i][j] + dU3dt_bar[i][j]) / 2;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				dU5dt_av[i][j] = (dU5dt[i][j] + dU5dt_bar[i][j]) / 2;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U1[i][j] += dU1dt_av[i][j] * dt;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U2[i][j] += dU2dt_av[i][j] * dt;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U3[i][j] += dU3dt_av[i][j] * dt;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U5[i][j] += dU5dt_av[i][j] * dt;
			}
		}
		getPrimitiveVariables(rho, u, v, T, p, e, mu, k, Ma, U1, U2, U3, U5);
		boundaryCondition(rho, u, v, T, p, e, mu, k, Ma, U1, U2, U3, U5);
	}
	void boundaryCondition(std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& e, std::vector<std::vector<double>>& mu, std::vector<std::vector<double>>& k, std::vector<std::vector<double>>& Ma, std::vector<std::vector<double>>& U1, std::vector<std::vector<double>>& U2, std::vector<std::vector<double>>& U3, std::vector<std::vector<double>>& U5) {
		u[0][0] = 0;
		v[0][0] = 0;
		p[0][0] = p_inf;
		T[0][0] = T_inf;
		for (int i = 1; i < Ny + 1; i++) {
			u[0][i] = Ma_inf * sqrt(gamma * R * T_inf);
			v[0][i] = 0;
			p[0][i] = p_inf;
			T[0][i] = T_inf;
		}
		for (int i = 1; i < Nx + 1; i++) {
			u[i][Ny] = Ma_inf * sqrt(gamma * R * T_inf);
			v[i][Ny] = 0;
			p[i][Ny] = p_inf;
			T[i][Ny] = T_inf;
		}
		for (int i = 1; i < Nx + 1; i++) {
			u[i][0] = 0;
			v[i][0] = 0;
			p[i][0] = 2 * p[i][1] - p[i][2];
			if (wall == "isothermal") {
				T[i][0] = T_w;
			}
			else if (wall == "adiabatic") {
				T[i][0] = (4 * T[i][1] - T[i][2]) / 3;
			}
		}
		for (int i = 1; i < Ny; i++) {
			u[Nx][i] = 2 * u[Nx - 1][i] - u[Nx - 2][i];
			v[Nx][i] = 2 * v[Nx - 1][i] - v[Nx - 2][i];
			p[Nx][i] = 2 * p[Nx - 1][i] - p[Nx - 2][i];
			T[Nx][i] = 2 * T[Nx - 1][i] - T[Nx - 2][i];
		}
		for (int i = 0; i < Ny + 1; i++) {
			rho[0][i] = p[0][i] / (R * T[0][i]);
			e[0][i] = c_v * T[0][i];
			mu[0][i] = mu0 * pow(T[0][i] / T0, 1.5) * (T0 + 110) / (T[0][i] + 110);
			k[0][i] = mu[0][i] * c_p / Pr;
			Ma[0][i] = sqrt(pow(u[0][i], 2) + pow(v[0][i], 2)) / sqrt(gamma * R * T[0][i]);
			U1[0][i] = rho[0][i];
			U2[0][i] = rho[0][i] * u[0][i];
			U3[0][i] = rho[0][i] * v[0][i];
			U5[0][i] = rho[0][i] * (e[0][i] + (pow(u[0][i], 2) + pow(v[0][i], 2)) / 2);
		}
		for (int i = 0; i < Ny + 1; i++) {
			rho[Nx][i] = p[Nx][i] / (R * T[Nx][i]);
			e[Nx][i] = c_v * T[Nx][i];
			mu[Nx][i] = mu0 * pow(T[Nx][i] / T0, 1.5) * (T0 + 110) / (T[Nx][i] + 110);
			k[Nx][i] = mu[Nx][i] * c_p / Pr;
			Ma[Nx][i] = sqrt(pow(u[Nx][i], 2) + pow(v[Nx][i], 2)) / sqrt(gamma * R * T[Nx][i]);
			U1[Nx][i] = rho[Nx][i];
			U2[Nx][i] = rho[Nx][i] * u[Nx][i];
			U3[Nx][i] = rho[Nx][i] * v[Nx][i];
			U5[Nx][i] = rho[Nx][i] * (e[Nx][i] + (pow(u[Nx][i], 2) + pow(v[Nx][i], 2)) / 2);
		}
		for (int i = 1; i < Nx; i++) {
			rho[i][0] = p[i][0] / (R * T[i][0]);
			e[i][0] = c_v * T[i][0];
			mu[i][0] = mu0 * pow(T[i][0] / T0, 1.5) * (T0 + 110) / (T[i][0] + 110);
			k[i][0] = mu[i][0] * c_p / Pr;
			Ma[i][0] = sqrt(pow(u[i][0], 2) + pow(v[i][0], 2)) / sqrt(gamma * R * T[i][0]);
			U1[i][0] = rho[i][0];
			U2[i][0] = rho[i][0] * u[i][0];
			U3[i][0] = rho[i][0] * v[i][0];
			U5[i][0] = rho[i][0] * (e[i][0] + (pow(u[i][0], 2) + pow(v[i][0], 2)) / 2);
		}
		for (int i = 1; i < Nx; i++) {
			rho[i][Ny] = p[i][Ny] / (R * T[i][Ny]);
			e[i][Ny] = c_v * T[i][Ny];
			mu[i][Ny] = mu0 * pow(T[i][Ny] / T0, 1.5) * (T0 + 110) / (T[i][Ny] + 110);
			k[i][Ny] = mu[i][Ny] * c_p / Pr;
			Ma[i][Ny] = sqrt(pow(u[i][Ny], 2) + pow(v[i][Ny], 2)) / sqrt(gamma * R * T[i][Ny]);
			U1[i][Ny] = rho[i][Ny];
			U2[i][Ny] = rho[i][Ny] * u[i][Ny];
			U3[i][Ny] = rho[i][Ny] * v[i][Ny];
			U5[i][Ny] = rho[i][Ny] * (e[i][Ny] + (pow(u[i][Ny], 2) + pow(v[i][Ny], 2)) / 2);
		}
	}
	double residual(Data& data) {
		double res = abs(data.rho[data.rho.size() - 1][0][0] - data.rho[data.rho.size() - 2][0][0]);
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				res = res < abs(data.rho[data.rho.size() - 1][i][j] - data.rho[data.rho.size() - 2][i][j]) ? abs(data.rho[data.rho.size() - 1][i][j] - data.rho[data.rho.size() - 2][i][j]) : res;
			}
		}
		return res;
	}
	void getU(std::vector<std::vector<double>>& U1, std::vector<std::vector<double>>& U2, std::vector<std::vector<double>>& U3, std::vector<std::vector<double>>& U5, std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& e) {
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U1[i][j] = rho[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U2[i][j] = rho[i][j] * u[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U3[i][j] = rho[i][j] * v[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				U5[i][j] = rho[i][j] * (e[i][j] + (pow(u[i][j], 2) + pow(v[i][j], 2)) / 2);
			}
		}
	}
	void getE(std::vector<std::vector<double>>& E1, std::vector<std::vector<double>>& E2, std::vector<std::vector<double>>& E3, std::vector<std::vector<double>>& E5, std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& e, std::vector<std::vector<double>>& tau_xx, std::vector<std::vector<double>>& tau_xy, std::vector<std::vector<double>>& q_x) {
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				E1[i][j] = rho[i][j] * u[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				E2[i][j] = rho[i][j] * pow(u[i][j], 2) + p[i][j] - tau_xx[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				E3[i][j] = rho[i][j] * u[i][j] * v[i][j] - tau_xy[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				E5[i][j] = (rho[i][j] * (e[i][j] + (pow(u[i][j], 2) + pow(v[i][j], 2)) / 2) + p[i][j]) * u[i][j] - u[i][j] * tau_xx[i][j] - v[i][j] * tau_xy[i][j] + q_x[i][j];
			}
		}
	}
	void getF(std::vector<std::vector<double>>& F1, std::vector<std::vector<double>>& F2, std::vector<std::vector<double>>& F3, std::vector<std::vector<double>>& F5, std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& e, std::vector<std::vector<double>>& tau_xy, std::vector<std::vector<double>>& tau_yy, std::vector<std::vector<double>>& q_y) {
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				F1[i][j] = rho[i][j] * v[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				F2[i][j] = rho[i][j] * u[i][j] * v[i][j] - tau_xy[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				F3[i][j] = rho[i][j] * pow(v[i][j], 2) + p[i][j] - tau_yy[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				F5[i][j] = (rho[i][j] * (e[i][j] + (pow(u[i][j], 2) + pow(v[i][j], 2)) / 2) + p[i][j]) * v[i][j] - u[i][j] * tau_xy[i][j] - v[i][j] * tau_yy[i][j] + q_y[i][j];
			}
		}
	}
	void getPrimitiveVariables(std::vector<std::vector<double>>& rho, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& p, std::vector<std::vector<double>>& e, std::vector<std::vector<double>>& mu, std::vector<std::vector<double>>& k, std::vector<std::vector<double>>& Ma, std::vector<std::vector<double>>& U1, std::vector<std::vector<double>>& U2, std::vector<std::vector<double>>& U3, std::vector<std::vector<double>>& U5) {
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				rho[i][j] = U1[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				u[i][j] = U2[i][j] / U1[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				v[i][j] = U3[i][j] / U1[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				e[i][j] = U5[i][j] / U1[i][j] - (pow(u[i][j], 2) + pow(v[i][j], 2)) / 2;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				T[i][j] = e[i][j] / c_v;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				p[i][j] = rho[i][j] * R * T[i][j];
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				mu[i][j] = mu0 * pow(T[i][j] / T0, 1.5) * (T0 + 110) / (T[i][j] + 110);
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				k[i][j] = mu[i][j] * c_p / Pr;
			}
		}
		for (int i = 0; i < Nx + 1; i++) {
			for (int j = 0; j < Ny + 1; j++) {
				Ma[i][j] = sqrt(pow(u[i][j], 2) + pow(v[i][j], 2)) / sqrt(gamma * R * T[i][j]);
			}
		}
	}
	void getTauxx(std::vector<std::vector<double>>& tau_xx, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& mu, int mode) {
		switch (mode) {
		case 0:
			for (int i = 1; i < Nx + 1; i++) {
				for (int j = 1; j < Ny; j++) {
					tau_xx[i][j] = -2 * mu[i][j] / 3 * ((u[i][j] - u[i - 1][j]) / dx + (v[i][j + 1] - v[i][j - 1]) / (2 * dy)) + 2 * mu[i][j] * (u[i][j] - u[i - 1][j]) / dx;
				}
			}
			break;
		case 1:
			for (int i = 0; i < Nx; i++) {
				for (int j = 1; j < Ny; j++) {
					tau_xx[i][j] = -2 * mu[i][j] / 3 * ((u[i + 1][j] - u[i][j]) / dx + (v[i][j + 1] - v[i][j - 1]) / (2 * dy)) + 2 * mu[i][j] * (u[i + 1][j] - u[i][j]) / dx;
				}
			}
			break;
		}
	}
	void getTauxy(std::vector<std::vector<double>>& tau_xy, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& mu, int mode) {
		switch (mode) {
		case 0:
			for (int i = 1; i < Nx + 1; i++) {
				for (int j = 1; j < Ny; j++) {
					tau_xy[i][j] = mu[i][j] * ((v[i][j] - v[i - 1][j]) / dx + (u[i][j + 1] - u[i][j - 1]) / (2 * dy));
				}
			}
			break;
		case 1:
			for (int i = 1; i < Nx; i++) {
				for (int j = 1; j < Ny + 1; j++) {
					tau_xy[i][j] = mu[i][j] * ((v[i + 1][j] - v[i - 1][j]) / (2 * dx) + (u[i][j] - u[i][j - 1]) / dy);
				}
			}
			break;
		case 2:
			for (int i = 0; i < Nx; i++) {
				for (int j = 1; j < Ny; j++) {
					tau_xy[i][j] = mu[i][j] * ((v[i + 1][j] - v[i][j]) / dx + (u[i][j + 1] - u[i][j - 1]) / (2 * dy));
				}
			}
			break;
		case 3:
			for (int i = 1; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					tau_xy[i][j] = mu[i][j] * ((v[i + 1][j] - v[i - 1][j]) / (2 * dx) + (u[i][j + 1] - u[i][j]) / dy);
				}
			}
			break;
		}
	}
	void getTauyy(std::vector<std::vector<double>>& tau_yy, std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, std::vector<std::vector<double>>& mu, int mode) {
		switch (mode) {
		case 0:
			for (int i = 1; i < Nx; i++) {
				for (int j = 1; j < Ny + 1; j++) {
					tau_yy[i][j] = -2 * mu[i][j] / 3 * ((u[i + 1][j] - u[i - 1][j]) / (2 * dx) + (v[i][j] - v[i][j - 1]) / dy) + 2 * mu[i][j] * (v[i][j] - v[i][j - 1]) / dy;
				}
			}
			break;
		case 1:
			for (int i = 1; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					tau_yy[i][j] = -2 * mu[i][j] / 3 * ((u[i + 1][j] - u[i - 1][j]) / (2 * dx) + (v[i][j + 1] - v[i][j]) / dy) + 2 * mu[i][j] * (v[i][j + 1] - v[i][j]) / dy;
				}
			}
			break;
		}
	}
	void getQx(std::vector<std::vector<double>>& q_x, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& k, int mode) {
		switch (mode) {
		case 0:
			for (int i = 1; i < Nx + 1; i++) {
				for (int j = 0; j < Ny + 1; j++) {
					q_x[i][j] = -k[i][j] * (T[i][j] - T[i - 1][j]) / dx;
				}
			}
			break;
		case 1:
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny + 1; j++) {
					q_x[i][j] = -k[i][j] * (T[i + 1][j] - T[i][j]) / dx;
				}
			}
			break;
		}
	}
	void getQy(std::vector<std::vector<double>>& q_y, std::vector<std::vector<double>>& T, std::vector<std::vector<double>>& k, int mode) {
		switch (mode) {
		case 0:
			for (int i = 0; i < Nx + 1; i++) {
				for (int j = 1; j < Ny + 1; j++) {
					q_y[i][j] = -k[i][j] * (T[i][j] - T[i][j - 1]) / dy;
				}
			}
			break;
		case 1:
			for (int i = 0; i < Nx + 1; i++) {
				for (int j = 0; j < Ny; j++) {
					q_y[i][j] = -k[i][j] * (T[i][j + 1] - T[i][j]) / dy;
				}
			}
			break;
		}
	}
public:
	SupersonicFlowOverAFlatPlate(double gamma, double R, double mu0, double T0, double Pr, double Ma_inf, double p_inf, double T_inf, double T_w, std::string wall, double LHORI, int Nx, int Ny, double K, double tol) :gamma(gamma), R(R), c_p(gamma* R / (gamma - 1)), c_v(R / (gamma - 1)), mu0(mu0), T0(T0), Pr(Pr), Ma_inf(Ma_inf), p_inf(p_inf), T_inf(T_inf), T_w(T_w), wall(wall), Nx(Nx), Ny(Ny), dx(LHORI / Nx), dy(25 * LHORI / (Ny * sqrt(p_inf / (R * T_inf) * Ma_inf * sqrt(gamma * R * T_inf) * LHORI / (mu0 * pow(T_inf / T0, 1.5) * (T0 + 110) / (T_inf + 110))))), K(K), tol(tol) {}
	void solve(Data& data) {
		double dt = 0;
		std::vector<std::vector<double>> rho, u, v, T, p, e, mu, k, Ma, U1, U2, U3, U5, tau_xx0, tau_xx1, tau_xy0, tau_xy1, tau_xy2, tau_xy3, tau_yy0, tau_yy1, q_x0, q_x1, q_y0, q_y1, E1, E2, E3, E5, F1, F2, F3, F5, dU1dt, dU2dt, dU3dt, dU5dt, U1_bar, U2_bar, U3_bar, U5_bar, rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, mu_bar, k_bar, Ma_bar, E1_bar, E2_bar, E3_bar, E5_bar, F1_bar, F2_bar, F3_bar, F5_bar, dU1dt_bar, dU2dt_bar, dU3dt_bar, dU5dt_bar, dU1dt_av, dU2dt_av, dU3dt_av, dU5dt_av;
		initialize(data, rho, u, v, T, p, e, mu, k, Ma, U1, U2, U3, U5, tau_xx0, tau_xx1, tau_xy0, tau_xy1, tau_xy2, tau_xy3, tau_yy0, tau_yy1, q_x0, q_x1, q_y0, q_y1, E1, E2, E3, E5, F1, F2, F3, F5, dU1dt, dU2dt, dU3dt, dU5dt, U1_bar, U2_bar, U3_bar, U5_bar, rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, mu_bar, k_bar, Ma_bar, E1_bar, E2_bar, E3_bar, E5_bar, F1_bar, F2_bar, F3_bar, F5_bar, dU1dt_bar, dU2dt_bar, dU3dt_bar, dU5dt_bar, dU1dt_av, dU2dt_av, dU3dt_av, dU5dt_av);
		save(data, dt, rho, u, v, T, p, Ma);
		getU(U1, U2, U3, U5, rho, u, v, e);
		while (true) {
			dt = timeStep(rho, u, v, T, mu);
			MacCormack(dt, rho, u, v, T, p, e, mu, k, Ma, U1, U2, U3, U5, tau_xx0, tau_xx1, tau_xy0, tau_xy1, tau_xy2, tau_xy3, tau_yy0, tau_yy1, q_x0, q_x1, q_y0, q_y1, E1, E2, E3, E5, F1, F2, F3, F5, dU1dt, dU2dt, dU3dt, dU5dt, U1_bar, U2_bar, U3_bar, U5_bar, rho_bar, u_bar, v_bar, T_bar, p_bar, e_bar, mu_bar, k_bar, Ma_bar, E1_bar, E2_bar, E3_bar, E5_bar, F1_bar, F2_bar, F3_bar, F5_bar, dU1dt_bar, dU2dt_bar, dU3dt_bar, dU5dt_bar, dU1dt_av, dU2dt_av, dU3dt_av, dU5dt_av);
			save(data, dt, rho, u, v, T, p, Ma);
			double res = residual(data);
			printf("time=%.8f res=%.8f\n", data.t[data.t.size() - 1], res);
			if (res < tol) {
				break;
			}
		}
	}
};
