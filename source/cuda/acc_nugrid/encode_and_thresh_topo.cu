#include "encode_and_thresh_topo.cuh"

__global__ void lis::cuda::acc_nugrid::encode_and_thresh_topo
(
	ScaleCoefficients d_scale_coeffs,
	Details           d_details,
	bool*             d_sig_details,
	Maxes             maxes,
	NUMERIC_TYPE      epsilon_local,
	int               level,
	bool              non_uniform_n,
	int               startfile
)
{
	index_1D idx = blockIdx.x * blockDim.x + threadIdx.x;

	int num_threads = 1 << (2 * level);

	if (idx >= num_threads) return;

	index_1D prev_lvl_idx = get_lvl_idx(level - 1);
	index_1D curr_lvl_idx = get_lvl_idx(level);
	index_1D next_lvl_idx = get_lvl_idx(level + 1);

	index_1D parent_idx = curr_lvl_idx + idx;
	index_1D child_idx  = next_lvl_idx + 4 * idx;

	ScaleChildren z0_children =
	{
		d_scale_coeffs.z0[child_idx + 0],
		d_scale_coeffs.z0[child_idx + 1],
		d_scale_coeffs.z0[child_idx + 2],
		d_scale_coeffs.z0[child_idx + 3]
	};

	ScaleChildren z1x_children =
	{
		d_scale_coeffs.z1x[child_idx + 0],
		d_scale_coeffs.z1x[child_idx + 1],
		d_scale_coeffs.z1x[child_idx + 2],
		d_scale_coeffs.z1x[child_idx + 3]
	};

	ScaleChildren z1y_children =
	{
		d_scale_coeffs.z1y[child_idx + 0],
		d_scale_coeffs.z1y[child_idx + 1],
		d_scale_coeffs.z1y[child_idx + 2],
		d_scale_coeffs.z1y[child_idx + 3]
	};

	SubDetail z_details =
	{
		encode_detail_alpha_11(z0_children,z1x_children,z1y_children),
		encode_detail_alpha_21(z0_children,z1x_children,z1y_children),
		encode_detail_alpha_12(z0_children,z1x_children,z1y_children),
		encode_detail_beta_11(z0_children,z1x_children,z1y_children),
		encode_detail_beta_21(z0_children,z1x_children,z1y_children),
		encode_detail_beta_12(z0_children,z1x_children,z1y_children),
		encode_detail_gamma_11(z0_children,z1x_children,z1y_children),
		encode_detail_gamma_21(z0_children,z1x_children,z1y_children),
		encode_detail_gamma_12(z0_children,z1x_children,z1y_children)
	};

	d_details.z0.alpha[parent_idx] = z_details.alpha_11;
	d_details.z0.beta[parent_idx]  = z_details.beta_11;
	d_details.z0.gamma[parent_idx] = z_details.gamma_11;

	d_details.z1x.alpha[parent_idx] = z_details.alpha_21;
	d_details.z1x.beta[parent_idx] = z_details.beta_21;
	d_details.z1x.gamma[parent_idx] = z_details.gamma_21;

	d_details.z1y.alpha[parent_idx] = z_details.alpha_12;
	d_details.z1y.beta[parent_idx] = z_details.beta_12;
	d_details.z1y.gamma[parent_idx] = z_details.gamma_12;

	d_scale_coeffs.z0[parent_idx] = encode_scale_11(z0_children, z1x_children, z1y_children);
	d_scale_coeffs.z1x[parent_idx] = encode_scale_21(z0_children, z1x_children, z1y_children);
	d_scale_coeffs.z1y[parent_idx] = encode_scale_12(z0_children, z1x_children, z1y_children);


	if (non_uniform_n) {
		ScaleChildren n0_children =
		{
			d_scale_coeffs.n0[child_idx + 0],
			d_scale_coeffs.n0[child_idx + 1],
			d_scale_coeffs.n0[child_idx + 2],
			d_scale_coeffs.n0[child_idx + 3]
		};

		ScaleChildren n1x_children =
		{
			C(0.0),
			C(0.0),
			C(0.0),
			C(0.0)
		};

		ScaleChildren n1y_children =
		{
			C(0.0),
			C(0.0),
			C(0.0),
			C(0.0)
		};

		SubDetail n_details =
		{
			encode_detail_alpha_11(n0_children,n1x_children,n1y_children),
			encode_detail_alpha_21(n0_children,n1x_children,n1y_children),
			encode_detail_alpha_12(n0_children,n1x_children,n1y_children),
			encode_detail_beta_11(n0_children,n1x_children,n1y_children),
			encode_detail_beta_21(n0_children,n1x_children,n1y_children),
			encode_detail_beta_12(n0_children,n1x_children,n1y_children),
			encode_detail_gamma_11(n0_children,n1x_children,n1y_children),
			encode_detail_gamma_21(n0_children,n1x_children,n1y_children),
			encode_detail_gamma_12(n0_children,n1x_children,n1y_children)
		};

		d_details.n0.alpha[parent_idx] = n_details.alpha_11;
		d_details.n0.beta[parent_idx] = n_details.beta_11;
		d_details.n0.gamma[parent_idx] = n_details.gamma_11;

		d_details.n1x.alpha[parent_idx] = n_details.alpha_21;
		d_details.n1x.beta[parent_idx] = n_details.beta_21;
		d_details.n1x.gamma[parent_idx] = n_details.gamma_21;

		d_details.n1y.alpha[parent_idx] = n_details.alpha_12;
		d_details.n1y.beta[parent_idx] = n_details.beta_12;
		d_details.n1y.gamma[parent_idx] = n_details.gamma_12;

		d_scale_coeffs.n0[parent_idx] = encode_scale_11(n0_children, n1x_children, n1y_children);
//		d_scale_coeffs.n1x[parent_idx] = encode_scale_21(z0_children, z1x_children, z1y_children);
//		d_scale_coeffs.n1y[parent_idx] = encode_scale_12(z0_children, z1x_children, z1y_children);
	}

	if (startfile) {
		ScaleChildren h0_children =
		{
			d_scale_coeffs.h[child_idx + 0],
			d_scale_coeffs.h[child_idx + 1],
			d_scale_coeffs.h[child_idx + 2],
			d_scale_coeffs.h[child_idx + 3]
		};

		ScaleChildren h1x_children =
		{
			C(0.0),
			C(0.0),
			C(0.0),
			C(0.0)
		};

		ScaleChildren h1y_children =
		{
			C(0.0),
			C(0.0),
			C(0.0),
			C(0.0)
		};

		SubDetail h_details =
		{
			encode_detail_alpha_11(h0_children,h1x_children,h1y_children),
			encode_detail_alpha_21(h0_children,h1x_children,h1y_children),
			encode_detail_alpha_12(h0_children,h1x_children,h1y_children),
			encode_detail_beta_11(h0_children,h1x_children,h1y_children),
			encode_detail_beta_21(h0_children,h1x_children,h1y_children),
			encode_detail_beta_12(h0_children,h1x_children,h1y_children),
			encode_detail_gamma_11(h0_children,h1x_children,h1y_children),
			encode_detail_gamma_21(h0_children,h1x_children,h1y_children),
			encode_detail_gamma_12(h0_children,h1x_children,h1y_children)
		};

		d_details.h0.alpha[parent_idx] = h_details.alpha_11;
		d_details.h0.beta[parent_idx] = h_details.beta_11;
		d_details.h0.gamma[parent_idx] = h_details.gamma_11;

		d_details.h1x.alpha[parent_idx] = h_details.alpha_21;
		d_details.h1x.beta[parent_idx] = h_details.beta_21;
		d_details.h1x.gamma[parent_idx] = h_details.gamma_21;

		d_details.h1y.alpha[parent_idx] = h_details.alpha_12;
		d_details.h1y.beta[parent_idx] = h_details.beta_12;
		d_details.h1y.gamma[parent_idx] = h_details.gamma_12;

		d_scale_coeffs.h[parent_idx] = encode_scale_11(h0_children, h1x_children, h1y_children);
		//		d_scale_coeffs.n1x[parent_idx] = encode_scale_21(z0_children, z1x_children, z1y_children);
		//		d_scale_coeffs.n1y[parent_idx] = encode_scale_12(z0_children, z1x_children, z1y_children);

	}

	if ((z_details.get_max() / maxes.z) >= epsilon_local) {
	//	d_preflagged_details[parent_idx] = SIGNIFICANT;
		d_sig_details[parent_idx] = SIGNIFICANT;
	}
}