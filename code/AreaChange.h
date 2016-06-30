#ifndef AreaChange_H
#define AreaChange_H


class AreaChange {
    
    public:
                                    AreaChange(int a, double p);  
                                    AreaChange(AreaChange& a);
                                   ~AreaChange(void);
        bool                        operator<(const AreaChange& a) const;
        int                         getArea(void) { return area; }
        double                      getPosition(void) { return position; }
        void                        print(void);
    
    private:
        int                         area;
        double                      position;
};

class comp_history {
    
    public:
        bool                        operator()(AreaChange* m1, AreaChange* m2) const { return (*m1 < *m2); }
};

#endif