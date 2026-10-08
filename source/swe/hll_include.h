if (H_neg <= DepthThresh && H_pos <= DepthThresh)
{
	H_flux = C(0.0);
	HU_flux = C(0.0);
	HV_flux = C(0.0);
	return;
}

NUMERIC_TYPE U_neg, V_neg;
if (H_neg <= DepthThresh)
{
	U_neg = C(0.0);
	V_neg = C(0.0);
}
else
{
	U_neg = HU_neg / H_neg;
	V_neg = HV_neg / H_neg;
}

NUMERIC_TYPE U_pos, V_pos;
if (H_pos <= DepthThresh)
{
	U_pos = C(0.0);
	V_pos = C(0.0);
}
else
{
	U_pos = HU_pos / H_pos;
	V_pos = HV_pos / H_pos;
}

NUMERIC_TYPE A_neg = SQRT(g * H_neg);
NUMERIC_TYPE A_pos = SQRT(g * H_pos);

NUMERIC_TYPE H_star = pow(C(0.5)*(A_neg+A_pos) + C(0.25)*(U_neg - U_pos),
		C(2.0))	/ g;
NUMERIC_TYPE U_star = C(0.5)*(U_neg + U_pos) + A_neg - A_pos;
NUMERIC_TYPE A_star = SQRT(g*H_star);

NUMERIC_TYPE S_neg;
if (H_neg <= DepthThresh)
{
	S_neg = U_pos - C(2.0)*A_pos;
}
else
{
	S_neg = min(U_neg - A_neg, U_star - A_star);
}

NUMERIC_TYPE S_pos;
if (H_pos <= DepthThresh)
{
	S_pos = U_neg + C(2.0)*A_neg;
}
else
{
	S_pos = max(U_pos + A_pos, U_star + A_star);
}

NUMERIC_TYPE H_flux_neg = HU_neg;
NUMERIC_TYPE HU_flux_neg = U_neg*HU_neg + C(0.5)*g*H_neg*H_neg;
NUMERIC_TYPE HV_flux_neg = HV_neg*U_neg;

NUMERIC_TYPE H_flux_pos = HU_pos;
NUMERIC_TYPE HU_flux_pos = U_pos*HU_pos + C(0.5)*g*H_pos*H_pos;
NUMERIC_TYPE HV_flux_pos = HV_pos*U_pos;

if (S_neg >= C(0.0))
{
	H_flux = H_flux_neg;
	HU_flux = HU_flux_neg;
	HV_flux = HV_flux_neg;
}
else if (S_neg < C(0.0) && S_pos >= C(0.0))
{
	NUMERIC_TYPE H_flux_mid =
		(S_pos*H_flux_neg - S_neg*H_flux_pos + S_neg*S_pos*(H_pos-H_neg))
		/
		(S_pos - S_neg);
	NUMERIC_TYPE HU_flux_mid =
		(S_pos*HU_flux_neg - S_neg*HU_flux_pos + S_neg*S_pos*(HU_pos-HU_neg))
		/
		(S_pos - S_neg);

	H_flux = H_flux_mid;
	HU_flux = HU_flux_mid;

	NUMERIC_TYPE S_mid =
		(S_neg*H_pos*(U_pos - S_pos) - S_pos*H_neg*(U_neg - S_neg))
		/
		(H_pos*(U_pos-S_pos) - H_neg*(U_neg - S_neg));

	if (S_neg < C(0.0) && S_mid >= C(0.0))
	{
		HV_flux = H_flux_mid * V_neg;
	}
	else
	{
		HV_flux = H_flux_mid * V_pos;
	}
}
else
{
	H_flux = H_flux_pos;
	HU_flux = HU_flux_pos;
	HV_flux = HV_flux_pos;
}

