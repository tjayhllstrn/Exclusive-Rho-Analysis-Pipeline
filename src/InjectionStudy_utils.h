#ifndef InjectionStudy_utils_h
#define InjectionStudy_utils_h

struct inj_event{
    public:
    double phi_h;
    double Pol;
    double eps;
    double t;

    double truepipparent_pid;
    double truepho1parentparent_pid;
    double truepho2parentparent_pid;
    double truepho1parent_id;
    double truepho2parent_id;
    double truepho1parent_pid;
    double truepho2parent_pid;


    inj_event(double phi, double Pol, double eps, double t, 
        double truepipparent_pid, double truepho1parentparent_pid, double truepho2parentparent_pid, double truepho1parent_id, double truepho2parent_id, double truepho1parent_pid, double truepho2parent_pid
        ) : phi_h(phi), Pol(Pol), eps(eps), t(t) ,truepipparent_pid(truepipparent_pid), truepho1parentparent_pid(truepho1parentparent_pid), truepho2parentparent_pid(truepho2parentparent_pid), truepho1parent_id(truepho1parent_id), truepho2parent_id(truepho2parent_id), truepho1parent_pid(truepho1parent_pid), truepho2parent_pid(truepho2parent_pid) {}
};

class cross_section {
    public:
    double A0; // asymmetry intercept
    double At; // linear coefficient in t
    bool use_linear_t;

    cross_section(const char* asymmetry_type);

    double Asymmetry(inj_event event); // event-dependent asymmetry amplitude
    double Asymmetry(double t); // t-dependent asymmetry amplitude
    double AssignmentRule(inj_event event); // takes in the relevant kinematics and outputs whether the event is negative or positive helicity
    double CalcCC(inj_event event, double hel); // calculates the cross section for a given event and helicity
};


#endif