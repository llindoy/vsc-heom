#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <random>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <algorithm>
#include <chrono>

#include "../heom_operator_sparse.hpp"

#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include <linalg/decompositions/singular_value_decomposition/singular_value_decomposition.hpp>
#include "../molecular_system.hpp"
#include "../rkf.hpp"

#include "../bath_decomposition.hpp"



template <typename T>
void debye_bath_pade(const T& Lambda, const T& wc, const T& beta, size_t K, heom_bath<T>& terms)
{
    std::vector<T> gamma(K+1);
    using complex_type = typename correlation_terms<T>::complex_type;

    auto Ck = [&Lambda, &wc, &beta](T w){return 4*Lambda*wc/beta/(w*w-wc*wc);};

    std::vector<complex_type> S(K+1);
    complex_type A;
    gamma[0] = wc;
    S[0] =  Lambda*wc/std::tan(beta*wc/2.0);
    A =   -Lambda*wc;

    if(K != 0)
    {
        size_t N = K;
        pade_decomposition<T> decomp(N, beta);

        //now we bind the 
        for(size_t n=0; n < N; ++n)
        {
            gamma[N-n] = decomp.nu(n);
            S[N-n] = Ck(gamma[N-n])*gamma[N-n]*decomp.eta(n);
        }
    }
    for(size_t k = 0; k < gamma.size(); ++k){std::cerr << gamma[k] << std::endl;}
    std::cerr << "Using pade" << std::endl;
    terms.add_term(correlation_terms<T>(gamma, S, A));
}

template <typename T>
void debye_bath(const T& Lambda, const T& wc, const T& beta, size_t K, heom_bath<T>& terms, bool use_nz_trunc = false)
{
    std::vector<T> gamma(K+1);
    using complex_type = typename correlation_terms<T>::complex_type;
    std::vector<complex_type> S(K+1);
    complex_type A;

    //work out the frequencies
    gamma[0] = wc;
    for(size_t k=1; k<K+1; ++k)
    {
        gamma[k] = 2*M_PI*k/beta;
    }

    auto Ck = [&Lambda, &wc, &beta](T w){return 4*Lambda*wc/beta/(w*w-wc*wc);};
    
    //now work out the symmetric and antisymmetric contributions to the coupling constants
    
    S[0] =  Lambda*wc/std::tan(beta*wc/2.0);
    for(size_t k=1; k<K+1; ++k)
    {
        S[k] = Ck(gamma[k])*gamma[k];
    }

    A =     -Lambda*wc;

    terms.add_term(correlation_terms<T>(gamma, S, A));

    if(!use_nz_trunc)
    {
        T delta = 2.0*Lambda/(beta*wc);
        for(size_t k=0; k < K; ++k)
        {
            delta -= std::real(S[k])/std::real(gamma[k]);
        }
        std::cerr << "Using matsubara" << std::endl;

        terms.bind_truncation(standard_truncation<T>(delta));
    }
    else
    {
        size_t npsum = 3*K;
        pade_decomposition<T> decomp(npsum, beta);

        nakajima_zwanzig_truncation<T> nz;

        std::vector<T> _gamma(K-1);
        std::vector<complex_type> _S(K-1);
        for(size_t i = 0; i < K-1; ++i)
        {
            _gamma[i] = gamma[i+1];
            _S[i] = S[i+1];
        }
    
        nz.set_mats(correlation_terms<T>(_gamma, _S));

        std::vector<T> _gammap(npsum);
        std::vector<complex_type> _Sp(npsum);
        //now we bind the 
        for(size_t n=0; n < npsum; ++n)
        {
            _gammap[n] = decomp.nu(n);
            _Sp[n] = Ck(_gammap[n])*_gammap[n]*decomp.eta(n);
        }
        nz.set_pade(correlation_terms<T>(_gammap, _Sp));
        terms.bind_truncation(nz);
    }
}

int main(int argc, char* argv[])
{
    using real_type = double;
    using complex_type = linalg::complex<real_type>;
    using backend_type = linalg::blas_backend;
    using size_type = typename backend_type::size_type;
    backend_type::initialise();

    if(argc < 2)
    {
        std::cerr << argv[0] << " <input filename>" << std::endl;
        return 1;
    }

    try
    {
        std::ifstream ifs(argv[1]);
        if(!ifs.is_open())
        {
            std::cerr << "Could not open input file." << std::endl;
            return 1;
        }
    
        rapidjson::IStreamWrapper isw{ifs};

        rapidjson::Document doc {};
        doc.ParseStream(isw);

        if(doc.HasParseError())
        {
            std::cerr << "Error: " << doc.GetParseError() << std::endl << "Offset: " << doc.GetErrorOffset() << std::endl;
            return 1;
        }

        real_type from_cmn1 = 1.0/219474.63;
        real_type from_fs =  41.341374575751;
        
        dfs_replace_space_and_capitals_in_key(doc, doc.GetAllocator());

        molecular_system<real_type> system(doc["system"]);

        real_type beta = rapidjson_loader<real_type>::load(doc["beta"]);

        real_type tmax = rapidjson_loader<real_type>::load(doc["tmax"])*from_fs;
        real_type dt = rapidjson_loader<real_type>::load(doc["dt"])*from_fs;

        real_type wc = rapidjson_loader<real_type>::load(doc["wc"])*from_cmn1;

        real_type Lambda = 0;
        if(doc.HasMember("lambda"))
        {
            Lambda = rapidjson_loader<real_type>::load(doc["lambda"]);
        }
        else
        {
            ASSERT(doc.HasMember("alpha") && doc.HasMember("wb"), "No bath coupling strength specfied.");
            real_type fric = rapidjson_loader<real_type>::load(doc["alpha"]);
            real_type wb = rapidjson_loader<real_type>::load(doc["wb"])*from_cmn1;
            Lambda = wc*wb*fric/2.0;  
        }

        size_type L = rapidjson_loader<size_type>::load(doc["l"]);
        size_type K;


        bool use_pade = rapidjson_loader<bool>::load(doc["usepade"]);   

        bool use_cutoff = false;
        if(doc.HasMember("usecutoff"))
        {
            use_cutoff = rapidjson_loader<bool>::load(doc["usecutoff"]);   
        }

        //if(use_cutoff)
        //{
        //    K = std::ceil(beta*L*wc/(2*M_PI));
        //}
        //else
        //{
            K = rapidjson_loader<size_type>::load(doc["k"]);
        //}

        //construct the basis 
        system.construct_basis();

        //now build the system Hamiltonian
        linalg::matrix<real_type> ham;
        system.hamiltonian_operator(ham);

        //now build the position operator
        linalg::matrix<real_type> R2;
        system.R2_operator(R2);

        //now we add in the Caldiera Leggett renormalisation term to the Hamiltonian
        ham += Lambda*R2;

        linalg::matrix<real_type> R;
        system.position_operator(R);

        linalg::matrix<complex_type> Scoup(R);

        linalg::matrix<complex_type> Hsys(ham);

        linalg::matrix<real_type> rside;
        system.side_operator(rside);
        linalg::matrix<complex_type> side(rside);

        linalg::matrix<real_type> rrhobeta2;
        system.rho0(rrhobeta2, beta/2.0, Lambda);
        linalg::matrix<complex_type> rhobeta2(rrhobeta2);

        linalg::matrix<complex_type> temp(rrhobeta2.size(0), rrhobeta2.size(1));
        linalg::matrix<complex_type> rho0(rrhobeta2.size(0), rrhobeta2.size(1));
        temp.fill_zeros();
        rho0.fill_zeros();

        temp = side*rhobeta2;
        rho0 = rhobeta2*temp;


        complex_type zr = 0;
        for(size_t i = 0; i < rho0.size(0); ++i)
        {
            zr += rho0(i, i);
        }
        rho0/=zr;

        heom_bath<real_type> term;
        if(use_pade)
        {
            debye_bath_pade(Lambda, wc, beta, K, term);
        }
        else
        {
            debye_bath(Lambda, wc, beta, K, term, false);
        }
        heom_bath_operator_sparse<real_type> hop;
        hop.add_bath(term);

        if(!use_cutoff)
        {
            hop.initialise(L);
        }
        else
        {
            if(doc.HasMember("maxfactor"))
            {
                real_type maxfactor = rapidjson_loader<real_type>::load(doc["maxfactor"]);   
                if(maxfactor < 1){maxfactor = 1;}
                hop.initialise_cutoff(L*wc, maxfactor);
            }
            else
            {
                hop.initialise_cutoff(L*wc);
            }
        }
        std::cerr << hop.nados() << std::endl;

        size_t nhilb = Hsys.shape(0);

        linalg::vector<complex_type> A(hop.nados()*nhilb*nhilb);    A.fill_zeros();
        auto ados = A.reinterpret_shape(hop.nados(), nhilb, nhilb);
        ados[0] = rho0;

        rkf<complex_type> rk(ados.size());
        rk.dt() = dt;
        if(doc.HasMember("tol"))
        {
            rk.eps() =  rapidjson_loader<real_type>::load(doc["tol"]);
        }

        auto sidev = side.reinterpret_shape(nhilb*nhilb);
        auto rhov = ados[0].reinterpret_shape(nhilb*nhilb);

        std::cout << std::setprecision(16);
        std::cerr << std::setprecision(16);
        //thermalise
        {
            hop.set_interactions(Hsys, Scoup);

            auto start = std::chrono::high_resolution_clock::now();
            
            size_type nsteps = static_cast<size_t>(tmax/dt)+1;

            complex_type val = linalg::dot_product(linalg::conj(sidev), rhov);
            std::cout << "Running." << std::endl;
            std::cout << "0.0 " << linalg::real(val) << " 0 " << std::endl;
            //for(size_t j = 0; j < nhilb; ++j){ std::cout << std::real(ados[0](j, j)) << " ";}    std::cout << std::endl;


            for(size_t i = 0; i < nsteps; ++i)
            {
                rk(A, [&hop](const linalg::vector<complex_type>& ai, linalg::vector<complex_type>& oi){hop.apply(ai, oi, complex_type(0, -1));}, dt);


                val = linalg::dot_product(linalg::conj(sidev), rhov);

                auto curr = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(curr - start);
                std::cout << (i+1)*dt/from_fs << " " << linalg::real(val) << " " << duration.count() << std::endl;
                //for(size_t j = 0; j < nhilb; ++j){ std::cout << std::real(ados[0](j, j)) << " ";}    std::cout << std::endl;
            }
        }
        return 0;
    }
    catch(const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        RAISE_EXCEPTION("Failed to perform HEOM calculation.");
    }
}

