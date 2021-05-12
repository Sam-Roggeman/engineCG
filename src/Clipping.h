//
// Created by User on 8/05/2021.
//

#ifndef ENGINE_CLIPPING_H
#define ENGINE_CLIPPING_H

void calculatepFrBckPs(const Vector3D* X, const Vector3D* Y, const Vector3D* X2, const Vector3D* Y2, double& p, double& p2,const double dval ){
    p = (dval - Y->z) / (X->z - Y->z);
    p2 = (dval - Y2->z) / (X2->z - Y2->z);
}
void calculatepRLPs(const Vector3D* X, const Vector3D* Y, const Vector3D* X2, const Vector3D* Y2, double& p, double& p2,const double dval, const double dNear ){
    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
    p2 = (Y2->x * dNear+Y2->z*dval)/((Y2->x-X2->x)*dNear+(Y2->z-X2->z)*dval);
}
void calculatepUpDownPs(const Vector3D* X, const Vector3D* Y, const Vector3D* X2, const Vector3D* Y2, double& p, double& p2,const double dval, const double dNear ){
    p = (Y->y * dNear+Y->z*dval)/((Y->y-X->y)*dNear+(Y->z-X->z)*dval);
    p2 = (Y2->y * dNear+Y2->z*dval)/((Y2->y-X2->y)*dNear+(Y2->z-X2->z)*dval);
}


void clipPane(std::vector<Vector3D> pane, int i, std::vector<std::vector<Vector3D>>& l) {
    std::vector<std::vector<Vector3D>> nw_l = {};
    double dval;
    double dval2;
    bool a_inside ;
    bool b_inside ;
    bool c_inside ;
    const Vector3D* A;
    const Vector3D* B;
    const Vector3D* C;
    const Vector3D* X;
    const Vector3D* Y;
    const Vector3D* X2;
    const Vector3D* Y2;
    Vector3D A2;
    Vector3D B2;
    Vector3D C2;
    Vector3D D2;
    bool near = i == 0;
    bool far = i == 1;
    bool rechts = i ==2;
    bool links = i ==3;
    bool boven = i ==4;
    bool onder = i ==5;
    double dNear;
    double p;
    double p2;
    for ( auto triangle_it = l.begin(); triangle_it < l.end();triangle_it++){
        int count = 0;
        A = &triangle_it->at(0);
        B = &triangle_it->at(1);
        C = &triangle_it->at(2);
        if (near || far) {
            dval = pane[0].z;
            if (near) {
                a_inside = A->z <= dval;
                b_inside = B->z <= dval;
                c_inside = C->z <= dval;
            } else {
                a_inside = A->z >= dval;
                b_inside = B->z >= dval;
                c_inside = C->z >= dval;
            }
        }else if (rechts||links) {
            dNear = -std::max(std::max(pane[0].z, pane[1].z), std::max(pane[2].z, pane[3].z));
            dval = pane[0].x * (-dNear) / pane[0].z;
            if (rechts) {
                dval2 = A->x * (-dNear) / A->z;
                a_inside = (dval2 <= dval);
                dval2 = B->x * (-dNear) / B->z;
                b_inside = (dval2 <= dval);
                dval2 = C->x * (-dNear) / C->z;
                c_inside = (dval2 <= dval);
            } else {
                dval2 = A->x * (-dNear) / A->z;
                a_inside = (dval2 >= dval);
                dval2 = B->x * (-dNear) / B->z;
                b_inside = (dval2 >= dval);
                dval2 = C->x * (-dNear) / C->z;
                c_inside = (dval2 >= dval);
            }
        } else if (boven||onder) {
            dNear = -std::max(std::max(pane[0].z, pane[1].z), std::max(pane[2].z, pane[3].z));
            dval = -pane[0].y * dNear / pane[0].z;
            if (boven) {
                a_inside = -A->y * dNear / A->z <= dval;
                b_inside = -B->y * dNear / B->z <= dval;
                c_inside = -C->y * dNear / C->z <= dval;
            } else {
                a_inside = -A->y * dNear / A->z >= dval;
                b_inside = -B->y * dNear / B->z >= dval;
                c_inside = -C->y * dNear / C->z >= dval;
            }
        }
        if (a_inside) count++;
        if (b_inside) count++;
        if (c_inside) count++;
        if (count ==0) {
            continue;
        }
        if (count ==3) {
            A2 = *A;
            B2 = *B;
            C2 = *C;
            nw_l.emplace_back(std::vector<Vector3D>({A2,B2,C2}));
            continue;
        }
        else if (count == 1) {
            if (!a_inside && !b_inside) {
                //AC
                X = C; Y = A;
                //BC
                X2 = B; Y2 = C;
                C2 = *C;
            } else if (!a_inside && !c_inside) {
                //AB
                X = A; Y = B;
                //CB
                X2 = B; Y2 = C;
                C2 = *B;
            } else if (!c_inside && !b_inside) {
                //AB
                X = A; Y = B;
                //AC
                X2 = A; Y2 = C;
                C2 = *A;
            }
        }else if (count == 2){
            if (!a_inside) {
                //AB
                X = A; Y = B;
                //AC
                X2 = A; Y2 = C;
                B2 = *B;
                C2 = *C;
            } else if (!c_inside) {
                //AC
                X = C; Y = A;
                //BC
                X2 = B; Y2 = C;
                B2 = *A;
                C2 = *B;
            } else if (!b_inside) {
                //AB
                X = A; Y = B;
                //BC
                X2 = B; Y2 = C;
                B2 = *A;
                C2 = *C;
            }
        }

        if (near||far) {
            calculatepFrBckPs(X, Y, X2, Y2, p, p2, dval);
        }
        else if (boven||onder){
            calculatepUpDownPs(X,Y,X2,Y2,p,p2,dval,dNear);
        }
        else if (rechts||links){
            calculatepRLPs(X,Y,X2,Y2,p,p2,dval,dNear);
        }
        if (count == 2) {
            A2 = (p * *X) + ((1.0 - p) * *Y);
            D2 = (p2 * *X2) + ((1.0 - p2) * *Y2);
            nw_l.emplace_back(std::vector<Vector3D>({D2, A2, C2}));
            nw_l.emplace_back(std::vector<Vector3D>({A2,B2,C2}));
        }
        else if (count == 1){
            A2 = (p * *X) + ((1 - p) * *Y);
            B2 = (p2 * *X2) + ((1 - p2) * *Y2);
            nw_l.emplace_back(std::vector<Vector3D>({A2,B2,C2}));
        }
    }
    l.clear();
    l = nw_l;
}

std::vector<std::vector<Vector3D>> clipVistrum(const Vector3D &A ,const Vector3D &B,const Vector3D &C,const Figuur &vistrum,     std::vector<std::vector<Vector3D>>& l) {
    l.emplace_back(std::vector<Vector3D>({A,B,C}));
    for (int i = 0; i <6; i++){
        clipPane(vistrum.vlakInPoints(i), i,l);
    }
    return l;
}
#endif //ENGINE_CLIPPING_H
