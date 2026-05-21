#include "InjectionStudy_utils.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <limits>
#include <cstdlib>
#include <TRandom.h>

cross_section::cross_section(const char* asymmetry_type){
    const char* const_prefix = "const:";
    const size_t const_prefix_len = std::strlen(const_prefix);
    if (asymmetry_type && std::strncmp(asymmetry_type, const_prefix, const_prefix_len) == 0) {
        const char* value_str = asymmetry_type + const_prefix_len;
        char* endptr = nullptr;
        double value = std::strtod(value_str, &endptr);
        if (endptr != value_str) {
            std::cout << "VALUE: " << value << std::endl;
            A0 = value;
            At = 0.0;
            use_linear_t = false;
            return;
        }
    }

    const char* allowed_types[6];
    allowed_types[0] = "A01";
    allowed_types[1] = "A10";
    allowed_types[2] = "A20";
    allowed_types[3] = "Alin_t";
    allowed_types[4] = "Am04";
    allowed_types[5] = "A04";

    if(strcmp(asymmetry_type, "A01") == 0){
        A0 = 0.01;
        At = 0.0;
        use_linear_t = false;
    } else if (strcmp(asymmetry_type, "A10") == 0){
        A0 = 0.1;
        At = 0.0;
        use_linear_t = false;
    } else if (strcmp(asymmetry_type, "Am04") == 0){
        A0 = -0.04;
        At = 0.0;
        use_linear_t = false;
    } else if (strcmp(asymmetry_type, "A04") == 0){
        A0 = 0.04;
        At = 0.0;
        use_linear_t = false;
    } else if (strcmp(asymmetry_type, "A20") == 0){
        A0 = 0.2;
        At = 0.0;
        use_linear_t = false;
    } else if (strcmp(asymmetry_type, "Alin_t") == 0){
        A0 = 0.0;
        At = -0.1; // linear coefficient in t
        use_linear_t = true;
    } else {
        std::cerr << "\033[33mUnknown asymmetry type: " << asymmetry_type << ". choose from const:VALUE or \033[0m";
        for(int i = 0; i < sizeof(allowed_types)/sizeof(allowed_types[0]); i++){
            std::cerr << allowed_types[i] << " ";
        }
        std::cerr << std::endl;
        A0 = 0.0;
        At = 0.0;
        use_linear_t = false;
    }
}

double cross_section::Asymmetry(inj_event event){
    if (use_linear_t) {
        return A0 + At * event.t;
    }
    return A0;
}

double cross_section::Asymmetry(double t){
    if (use_linear_t) {
        return A0 + At * t;
    }
    return A0;
}

//Takes in the relevant kinematics and outputs whether the event is negative or positive helicity
double cross_section::AssignmentRule(inj_event event){
    double AssignedHel;
    const double sigma_pos = CalcCC(event, 1);
    const double sigma_neg = CalcCC(event, -1);

    if (sigma_neg <= std::numeric_limits<double>::epsilon()) {
        return 1.0;
    }

    double R = sigma_pos / sigma_neg; // calculate the ratio of cross sections for + vs - helicity
    double rand_num = gRandom ? gRandom->Uniform(0.0, 1.0) : 0.5;
    
    if(rand_num < R/(1+R)){
        AssignedHel = 1;
    } else {
        AssignedHel = -1;
    }
    return AssignedHel;
}

//This will calculate the corss section for a given eventy and asymmetry amplitude. This is based on a SIDIS process with the unpolarized modulations set to 0
double cross_section::CalcCC(inj_event event,double hel){
    double phi = event.phi_h;
    double Pol = event.Pol;
    double eps = event.eps;
    double depol = sqrt(2*eps*(1-eps));
    double Aevt = Asymmetry(event);

    double cc = 1 + hel*Pol*depol*Aevt*sin(phi);
    return cc;
}