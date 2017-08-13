#ifndef FASTMARCHING_HPP
#define FASTMARCHING_HPP

#include <queue>
#include <vector>
#include <cfloat>
#include <cmath>
#include <string>

#include "Debug.hpp"
#include "UtilityFunctions.hpp"

#define INF FLT_MAX
#define edge std::pair<unsigned int, float>

struct edgeComp {
    bool operator() (const edge &a, const edge &b) {
        return a.second > b.second;
    }
};

float FastMarching2dUpdate2P(float Ua, float Ub, float f)
{
    if (std::abs(Ua-Ub) <= 1.0/f)
    {
        return (Ua + Ub + std::sqrt(2*(1.0/std::pow(f, 2)) - std::pow(Ua - Ub, 2))) / 2.0;
    }
    else
    {
        return std::min(Ua, Ub) + 1.0/f;
    }
}

float FastMarching2dUpdate3P(float Ui, float Uj, float Uk, float f)
{
    if (Ui != INF and Uj != INF)
    {
        return FastMarching2dUpdate2P(Ui, Uj, f);
    }
    else if (Ui != INF and Uk != INF)
    {
        return FastMarching2dUpdate2P(Ui, Uk, f);
    }
    else if (Uj != INF and Uk != INF)
    {
        return FastMarching2dUpdate2P(Uj, Uk, f);
    }
    else
    {
        return std::min(Ui, std::min(Uj, Uk)) + 1.0/f;
    }
}

void FastMarching(int width, int height, int depth, float *sdm)
{
    Debug::Info("FastMarching: Entering");

    // Init values sizes and edges

    int size = width * height * depth;
    int numAcceptedNodes = 0;

    std::priority_queue< edge, std::vector<edge>, edgeComp> Q;

    int NidxX[26];
    int NidxY[26];
    int NidxZ[26];

    int idx = 0;
    for (int n = 1; n <= 3; ++n)
    {
        for (int k = -1; k <= 1; ++k)
        for (int j = -1; j <= 1; ++j)
        for (int i = -1; i <= 1; ++i)
        {
            if (std::abs(i) + std::abs(j) + std::abs(k) == n)
            {
                NidxX[idx] = i;
                NidxY[idx] = j;
                NidxZ[idx] = k;

                ++idx;
            }

        }
    }

    float *D = (float *) malloc(size * sizeof(float));

    if (D == NULL) {
        Debug::Error("FastMarching: Malloc failed (D)");
    }

    bool *F = (bool *) malloc(size * sizeof(bool));

    if (F == NULL) {
        Debug::Error("FastMarching: Malloc failed (F)");
    }

    for(int i = 0; i < size; i++)
    {
        if (sdm[i] == -2)
        {
            D[i] = 0;
            F[i] = false;
            Q.push(edge(i, 0));
        }
        else if (sdm[i] == -1)
        {
            D[i] = INF;
            F[i] = true;
        }
        else
        {
            D[i] = INF;
            F[i] = false;
        }
    }

    // Run Dijkstras Fast Marching

    Debug::Info("FastMarching: Running main loop");

    while(!Q.empty())
    {
        int u = Q.top().first;
        Q.pop();

        // already found
        if(F[u]) continue;

        numAcceptedNodes++;

        // mark as found
        F[u] = true;

        // DEBUG BEGIN

        int ui = u % width;
        int uj = (u / width) % height;
        int uk = u / (width*height);

        // DEBUG END

        // relax edges
        //int n = GN[u];
        //for(int i = 0; i < n; i++)
        for (int nindex = 0; nindex < 26; ++nindex)
        {
            int ni = NidxX[nindex];
            int nj = NidxY[nindex];
            int nk = NidxZ[nindex];

            int ii = ni + ui;
            int jj = nj + uj;
            int kk = nk + uk;

            if (ii < 0 or
                jj < 0 or
                kk < 0 or
                ii == width or
                jj == height or
                kk == depth)
            {
                continue;
            }


            int v = gIndex(ii, jj, kk, height, width, 0);

            if(!F[v])
            {
                float f = 1.0;

                bool threeDf = false;

                if (sdm[v] != -1 and ((jj == uj and kk == uk) or threeDf))
                {
                    if (sdm[v] == 0 or sdm[v] == -2)
                    {
                        f = 10000;
                    }
                    else
                    {
                        f = 1.0/sdm[v];
                    }

                }

                // Normal Dijkstra here
                //if (D[v] > D[u] + w[nindex]) D[v] = D[u] + w[nindex];

                if (std::abs(ni) + std::abs(nj) + std::abs(nk) > 1) continue;

                float Dip = INF, Dim = INF;
                float Djp = INF, Djm = INF;
                float Dkp = INF, Dkm = INF;

                if (ii > 0 and sdm[gIndex(ii-1,jj,kk,height,width,0)] != -1)          Dim = D[gIndex(ii-1,jj,kk,height,width,0)];
                if (ii < width - 1 and sdm[gIndex(ii+1,jj,kk,height,width,0)] != -1)  Dip = D[gIndex(ii+1,jj,kk,height,width,0)];
                if (jj > 0 and sdm[gIndex(ii,jj-1,kk,height,width,0)] != -1)          Djm = D[gIndex(ii,jj-1,kk,height,width,0)];
                if (jj < height - 1 and sdm[gIndex(ii,jj+1,kk,height,width,0)] != -1) Djp = D[gIndex(ii,jj+1,kk,height,width,0)];
                if (kk > 0 and sdm[gIndex(ii,jj,kk-1,height,width,0)] != -1)          Dkm = D[gIndex(ii,jj,kk-1,height,width,0)];
                if (kk < depth - 1 and sdm[gIndex(ii,jj,kk+1,height,width,0)] != -1)  Dkp = D[gIndex(ii,jj,kk+1,height,width,0)];

                float Ui = std::min(Dim, Dip);
                float Uj = std::min(Djm, Djp);
                float Uk = std::min(Dkm, Dkp);

                float newD;

                if (Ui == INF or Uj == INF or Uk == INF)
                {
                    if (Ui == INF and Uj == INF and Uk == INF)
                    {
                        // This happens if the node is completely isolated
                        // from the source points.
                        newD = INF;
                    }
                    else
                    {
                        newD = FastMarching2dUpdate3P(Ui, Uj, Uk, f);
                    }
                }
                else
                {
                    float sumU = Ui + Uj + Uk;
                    float d = sumU * sumU - 3*(Ui*Ui + Uj*Uj + Uk*Uk - 1.0/std::pow(f, 2));

                    if (d < 0)
                    {
                        newD = FastMarching2dUpdate3P(Ui, Uj, Uk, f);

                    }
                    else
                    {
                        newD = (sumU + std::sqrt(d)) / 3.0;
                    }
                }

                if (newD <= D[v])
                {
                    D[v] = newD;
                }
                Q.push(edge(v, D[v]));
            }

        }
    }

    for (int i = 0; i < size; ++i)
    {
        if (D[i] < INF and sdm[i] != -2)
        {
                sdm[i] = D[i];
        }
        else
        {
           sdm[i] = 0;
        }
    }

    free(D);
    free(F);

    Debug::Info("FastMarching: Leaving");
}

#endif // FASTMARCHING_HPP