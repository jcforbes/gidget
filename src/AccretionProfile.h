class FixedMesh;
class Debug;


class AccretionProfile {
    public:
        AccretionProfile(FixedMesh& mesh, int whichProfile, double alpha, Debug& dbg, double rs, double width);
        std::vector<double> & GetProfile() { return profile; };
        void UpdateProfile(double MhOverMh0);
        double fInner();
        double fOuter();
    private:
        FixedMesh & mesh;
        Debug & dbg;
        unsigned int nx;
        double radialScale, currentRadialScale;
        double alpha;
        int whichProfile;
        std::vector<double>  profile;
        double sigma, width;
        double normalization;

};



