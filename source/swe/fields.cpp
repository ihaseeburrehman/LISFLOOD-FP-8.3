#include "fields.h"
#include "../utility.h"

void allocate_swe_fields
(
	Pars *Parptr,
	Arrays *Arrptr
)
{
	Arrptr->HU = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HV = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->FHx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->FHUx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->FHVx = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->FHy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->FHUy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->FHVy = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->Zstar_x = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));
	Arrptr->Zstar_y = memory_allocate_zero_numeric_legacy((Parptr->xsz + 1)*(Parptr->ysz + 1));	
	Arrptr->Hstar_neg_x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->Hstar_pos_x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->Hstar_neg_y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->Hstar_pos_y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
}

void deallocate_swe_fields(Arrays *Arrptr)
{
	memory_free_legacy(&Arrptr->HU);
	memory_free_legacy(&Arrptr->HV);
	memory_free_legacy(&Arrptr->FHx);
	memory_free_legacy(&Arrptr->FHUx);
	memory_free_legacy(&Arrptr->FHVx);
	memory_free_legacy(&Arrptr->FHy);
	memory_free_legacy(&Arrptr->FHUy);
	memory_free_legacy(&Arrptr->FHVy);
	memory_free_legacy(&Arrptr->Zstar_x);
	memory_free_legacy(&Arrptr->Zstar_y);
	memory_free_legacy(&Arrptr->Hstar_neg_x);
	memory_free_legacy(&Arrptr->Hstar_pos_x);
	memory_free_legacy(&Arrptr->Hstar_neg_y);
	memory_free_legacy(&Arrptr->Hstar_pos_y);
}

