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

#include "../heom_operator.hpp"

#include <linalg/decompositions/eigensolvers/eigensolver.hpp>
#include <linalg/decompositions/singular_value_decomposition/singular_value_decomposition.hpp>
#include "../rkf.hpp"
#include "../molecular_system_2D.hpp"

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

        molecular_system_2D<real_type> system(doc["system"]);

        real_type beta = rapidjson_loader<real_type>::load(doc["beta"]);

        real_type tmax = rapidjson_loader<real_type>::load(doc["tmax"])*from_fs;
        real_type dt = rapidjson_loader<real_type>::load(doc["dt"])*from_fs;

        real_type wc = rapidjson_loader<real_type>::load(doc["wc"])*from_cmn1;

        real_type gammam = rapidjson_loader<real_type>::load(doc["gammam"])*from_cmn1;
        real_type Lambdam = 0;
        {
            real_type mode_lifetime = rapidjson_loader<real_type>::load(doc["modelifetime"])*from_fs;
            real_type loss_rate = 1.0/mode_lifetime;
            real_type w0 = system.mode_frequency();
            Lambdam = loss_rate/4.0 * (1-exp(-beta*w0)) * (gammam*gammam + w0*w0)/(w0*gammam)*w0;
        }

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
        size_type K = rapidjson_loader<size_type>::load(doc["k"]);
        size_type Km = rapidjson_loader<size_type>::load(doc["km"]);


        bool use_pade = rapidjson_loader<bool>::load(doc["usepade"]);   

        bool use_cutoff = false;
        if(doc.HasMember("usecutoff"))
        {
            use_cutoff = rapidjson_loader<bool>::load(doc["usecutoff"]);   
        }

        //construct the basis 
        system.construct_basis();

        //now build the system Hamiltonian
        linalg::matrix<real_type> ham;
        system.hamiltonian_operator(ham);

        //now build the position operator

        linalg::matrix<real_type> R2;
        system.R2_operator(R2);

        linalg::matrix<real_type> q2;
        system.Q2_operator(q2);

        //now we add in the Caldiera Leggett renormalisation term to the Hamiltonian

        linalg::matrix<real_type> R;
        system.R_operator(R);

        linalg::matrix<complex_type> Scoup(R);

        linalg::matrix<real_type> qc;
        system.Q_operator(qc);
        linalg::matrix<complex_type> Scoupc(qc);



        linalg::matrix<complex_type> Hsys0(ham);

        //also add on the bath renormalisation term
        ham += Lambda*R2;

        //and the molecule bath renormalisation term
        ham += Lambdam*q2;

        linalg::matrix<complex_type> Hsys(ham);

        linalg::matrix<real_type> rrhobeta2;
        system.rho0(rrhobeta2, beta/2.0, Lambda);
        linalg::matrix<complex_type> rhobeta2(rrhobeta2);

        linalg::matrix<complex_type> mu;

        bool dipole_R  = true;
        if(doc.HasMember("dipoler"))
        {
            dipole_R=  rapidjson_loader<bool>::load(doc["dipoler"]);
        }

        if(dipole_R)
        {
            mu = Scoup;
        }
        else
        {
            mu = Scoupc;
        }


        complex_type zr = 0;
        for(size_t i = 0; i < rhobeta2.size(0); ++i)
        {
            zr += rhobeta2(i, i);
        }
        rhobeta2/=zr;

        linalg::matrix<complex_type> rho0 = mu*rhobeta2;


        std::vector<heom_bath<real_type>> terms(2);
        if(use_pade)
        {
            debye_bath_pade(Lambda, wc, beta, K, terms[0]);
            debye_bath_pade(Lambdam, gammam, beta, Km, terms[1]);
        }
        else
        {
            debye_bath(Lambda, wc, beta, K, terms[0], false);
            debye_bath(Lambdam, gammam, beta, Km, terms[1], false);
        }
        heom_bath_operator<real_type> hop;
        hop.add_bath(terms[0]);
        hop.add_bath(terms[1]);


        if(!use_cutoff)
        {
            hop.initialise(L);
        }
        else
        {
            if(doc.HasMember("maxfactor"))
            {
                if(doc["maxfactor"].IsNumber())
                {
                    real_type maxfactor = rapidjson_loader<real_type>::load(doc["maxfactor"]);   
                    if(maxfactor < 1){maxfactor = 1;}
                    hop.initialise_cutoff(L*wc, maxfactor);
                }
                else
                {
                    if(doc["maxfactor"].IsArray())
                    {
                        ASSERT(doc["maxfactor"].Size() == 2, "Invalid max factor array");
                        for(size_t i = 0; i < 2; ++i){ASSERT(doc["maxfactor"][i].IsNumber(), "Invalid max factor array.");}
                        std::vector<real_type> maxfactor(2);
                        maxfactor[0] = rapidjson_loader<real_type>::load(doc["maxfactor"][0]);
                        maxfactor[1] = rapidjson_loader<real_type>::load(doc["maxfactor"][1]);

                        std::vector<real_type> scalefactor(2);  scalefactor[0] = 1; scalefactor[1] = 1;
                        try
                        {
                            if(doc.HasMember("scalefactor"))
                            {
                                ASSERT(doc["scalefactor"].IsArray(), "Invalid scalefactor.");
                                ASSERT(doc["scalefactor"].Size() == 2, "Invalid scalefactor.");
                                for(size_t i = 0; i < 2; ++i){ASSERT(doc["scalefactor"][i].IsNumber(), "Invalid scale factor array.");}
                                scalefactor[0] = rapidjson_loader<real_type>::load(doc["scalefactor"][0]);
                                scalefactor[1] = rapidjson_loader<real_type>::load(doc["scalefactor"][1]);
                            }
                        }
                        catch(const std::exception& ex)
                        {
                            std::cerr << "Failed to read scale factor object: " << std::endl << ex.what() << std::endl;
                        }
                        hop.initialise_cutoff(L*wc, maxfactor, scalefactor);
                    }
                }
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


        auto rhov = ados[0].reinterpret_shape(nhilb*nhilb);
        auto Rv = mu.reinterpret_shape(nhilb*nhilb);

        std::cout << std::setprecision(16);
        std::cerr << std::setprecision(16);
        
        //thermalise
        {
            auto Hop = superoperator::factory::construct(Hsys);
            std::vector<std::shared_ptr<superoperator::superoperator<complex_type>>> Sops(2);
            Sops[0] = superoperator::factory::construct(Scoup);
            Sops[1] = superoperator::factory::construct(Scoupc);
            
            linalg::csr_matrix<complex_type> M;
            hop(Hop, Sops, M);
            M.prune();

            size_type nsteps = static_cast<size_t>(tmax/dt)+1;

            auto start = std::chrono::high_resolution_clock::now();
            complex_type val0 = linalg::dot_product(linalg::conj(Rv), rhov);
            std::cout << "Running." << std::endl;
            std::cout << "0.0 " << linalg::real(val0/std::abs(val0)) << " " << linalg::imag(val0/std::abs(val0)) << " " << 0 << std::endl;
            //for(size_t j = 0; j < nhilb; ++j){ std::cout << std::real(ados[0](j, j)) << " ";}    std::cout << std::endl;

            linalg::diagonal_matrix<real_type> S;
            linalg::singular_value_decomposition<linalg::matrix<complex_type> > svd;


            for(size_t i = 0; i < nsteps; ++i)
            {
                rk(A, [&M](const linalg::vector<complex_type>& ai, linalg::vector<complex_type>& oi){oi = complex_type(0, -1)*M*ai;}, dt);

                complex_type val = linalg::dot_product(linalg::conj(Rv), rhov);

                auto curr = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(curr - start);
                std::cout << (i+1)*dt/from_fs << " " << linalg::real(val/std::abs(val0)) << " " << linalg::imag(val/std::abs(val0)) << " " << duration.count() << std::endl;
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


