#ifndef CShore_H
#define CShore_H
extern "C"
{
    void CShore(int *);

    void CShoreWrapper(int const*,        /* In_ILINE */
                       int const*,        /* In_IPROFL */
                       int const*,        /* In_IPERM */
                       int const*,        /* In_IOVER */
                       int const*,        /* In_IWCINT */
                       int const*,        /* In_IROLL */
                       int const*,        /* In_IWIND */
                       int const*,        /* In_ITIDE */
                       int const*,        /* In_ILAB */
                       int const*,        /* In_NWAVE */
                       int const*,        /* In_NSURG */
                       double const*,     /* In_DX */
                       double const*,     /* In_GAMMA */
                       double[],          /* In_TWAVE */
                       double[],          /* In_TPIN */
                       double[],          /* In_HRMSIN */
                       double[],          /* In_WANGIN */
                       double[],          /* In_TSURG */
                       double[],          /* In_SWLIN */
                       int const*,        /* In_NBINP */
                       double[],          /* In_XBINP */
                       double[],          /* In_ZBINP */
                       double[],          /* In_FBINP */
                       int*,              /* Out_IError */
                       int*,              /* Out_nOutSize */
                       double[],          /* Out_XYDist */
                       double[],          /* Out_FreeSurfaceStd */
                       double[],          /* Out_WaveSetupSurge */
                       double[],          /* Out_SinWaveAngleRadians */
                       double[]);         /* Out_FractionBreakingWaves */
}
#endif // CShore_H
