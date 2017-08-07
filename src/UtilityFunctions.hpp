#ifndef UTILITYFUNCTIONS_HPP
#define UTILITYFUNCTIONS_HPP

//#define gIndex(x,y,z,ny,nx,s) (z+s)*(ny+2*s)*(nx+2*s) + (y+s)*(nx+2*s) + x+s

#define STR(x) std::to_string(x)
#include <sstream> //for std::stringstream
#include <string>  //for std::string

inline int gIndex(int x, int y, int nx)
{
    return y*nx + x;
}

inline int gIndex(int x, int y, int nx, int s)
{
    return (y+s)*(nx+2*s) + x+s;
}

inline int gIndex(int x, int y, int z, int ny, int nx)
{
    return z*ny*nx + y*nx + x;
}

inline int gIndex(int x, int y, int z, int ny, int nx, int s)
{
    return (z+s)*(ny+2*s)*(nx+2*s) + (y+s)*(nx+2*s) + x+s;
}

inline int gIndex(int x, int y, int z, int k, int nz, int ny, int nx)
{
    return k*nz*ny*nx + z*ny*nx + y*nx + x;
}

template <class PTR_t>
std::string ASTR(PTR_t address)
{
    const void * vaddress = static_cast<const void*>(address);
    std::stringstream ss;
    ss << vaddress;
    return ss.str();
}

template <class T1, class T2, class T3>
T1 clamp(T1 v, T2 a, T3 b)
{
    return std::max(std::min(v, (T1) b), (T1) a);
}

#endif // UTILITYFUNCTIONS_HPP