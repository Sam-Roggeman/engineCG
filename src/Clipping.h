//
// Created by User on 8/05/2021.
//

#ifndef ENGINE_CLIPPING_H
#define ENGINE_CLIPPING_H
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
            if (a_inside) count++;
            if (b_inside) count++;
            if (c_inside) count++;
            if (count ==3) {
                A2 = *A;
                B2 = *B;
                C2 = *C;
            } else if (count ==0) {
                break;
            } else if (count == 1) {
                if (!a_inside && !b_inside) {
                    //AC
                    p = (dval - C->z) / (A->z - C->z);
                    A2 = (p * *A) + ((1 - p) * *C);
                    //BC
                    p = (dval - C->z) / (B->z - C->z);
                    B2 = (p * *B) + ((1 - p) * *C);
                    C2 = *C;
                } else if (!a_inside && !c_inside) {
                    //AB
                    p = (dval - B->z) / (A->z - B->z);
                    A2 = (p * *A) + ((1 - p) * *B);
                    //CB
                    p = (dval - B->z) / (C->z - B->z);
                    B2 = (p * *C) + ((1 - p) * *B);
                    C2 = *B;
                } else if (!c_inside && !b_inside) {
                    //CA
                    p = (dval - C->z) / (A->z - C->z);
                    A2 = (p * *C) + ((1 - p) * *A);
                    //BA
                    p = (dval - A->z) / (B->z - A->z);
                    B2 = (p * *B) + ((1 - p) * *A);
                    C2 = *A;
                }
            }else if (count == 2){
                if (!a_inside) {
                    //AC
                    X = A; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AB
                    X = A; Y = B;
                    p = (dval - Y->z) / (X->z - Y->z);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *C;
                    C2 = *B;
                } else if (!c_inside) {
                    //AC
                    X = A; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *B;
                } else if (!b_inside) {
                    //AB
                    X = A; Y = B;
                    p = (dval - Y->z) / (X->z - Y->z);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (dval - Y->z) / (X->z - Y->z);
                    D2 = (p * *X) + ((1 - p) * *Y);

                    B2 = *A;
                    C2 = *C;
                }
            }
        }else if (rechts||links){
            dNear = -std::max(std::max(pane[0].z, pane[1].z),std::max(pane[2].z, pane[3].z));
            dval = pane[0].x*(-dNear)/pane[0].z;
            if (rechts) {
                dval2 = A->x*(-dNear)/A->z;
                a_inside = (dval2 <= dval);
                dval2 = B->x*(-dNear)/B->z;
                b_inside = (dval2 <= dval);
                dval2 = C->x*(-dNear)/C->z;
                c_inside = (dval2 <= dval);
            } else {
                dval2 = A->x*(-dNear)/A->z;
                a_inside = (dval2 >= dval);
                dval2 = B->x*(-dNear)/B->z;
                b_inside = (dval2 >= dval);
                dval2 = C->x*(-dNear)/C->z;
                c_inside = (dval2 >= dval);
            }
            if (a_inside) count++;
            if (b_inside) count++;
            if (c_inside) count++;
            if (count ==3) {
                A2 = *A;
                B2 = *B;
                C2 = *C;
            } else if (count ==0) {
                break;
            }else if (count == 1) {
                if (!a_inside && !b_inside) {
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *C;
                } else if (!a_inside && !c_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //CB
                    X = C; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *B;
                } else if (!c_inside && !b_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *A;
                }
            }else if (count == 2){
                if (!a_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *B;
                    C2 = *C;
                } else if (!c_inside) {
                    //AC
                    X = A; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *B;
                } else if (!b_inside) {
                    //AB
                    X = A; Y = B;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B; Y = C;
                    p = (Y->x * dNear+Y->z*dval)/((Y->x-X->x)*dNear+(Y->z-X->z)*dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *C;
                }
            }
        } else if (boven||onder) {
            dNear = -std::max(std::max(pane[0].z, pane[1].z), std::max(pane[2].z, pane[3].z));
            dval = -pane[0].y * dNear / pane[0].z;
            if (boven) {
                a_inside =  -A->y*dNear/A->z <= dval;
                b_inside = -B->y*dNear/B->z <= dval;
                c_inside = -C->y*dNear/C->z <= dval;
            } else {
                a_inside = -A->y*dNear/A->z >= dval;
                b_inside = -B->y*dNear/B->z >= dval;
                c_inside = -C->y*dNear/C->z >= dval;
            }
            if (a_inside) count++;
            if (b_inside) count++;
            if (c_inside) count++;
            if (count ==3) {
                A2 = *A;
                B2 = *B;
                C2 = *C;
            } else if (count ==0) {
                break;
            } else if (count == 1) {
                if (!a_inside && !b_inside) {
                    //AC
                    X = A;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *C;
                } else if (!a_inside && !c_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //CB
                    X = C;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *B;
                } else if (!c_inside && !b_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    B2 = (p * *X) + ((1 - p) * *Y);
                    C2 = *A;
                }
            } else if (count == 2) {
                if (!a_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //AC
                    X = A;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *B;
                    C2 = *C;
                } else if (!c_inside) {
                    //BC
                    X = B;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //CA
                    X = C;
                    Y = A;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *B;
                    C2 = *A;
                } else if (!b_inside) {
                    //AB
                    X = A;
                    Y = B;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    A2 = (p * *X) + ((1 - p) * *Y);
                    //BC
                    X = B;
                    Y = C;
                    p = (Y->y * dNear + Y->z * dval) / ((Y->y - X->y) * dNear + (Y->z - X->z) * dval);
                    D2 = (p * *X) + ((1 - p) * *Y);
                    B2 = *A;
                    C2 = *C;
                }
            }
        }
        if (count == 2) {
            nw_l.emplace_back(std::vector<Vector3D>({D2, A2, C2}));
        }
        nw_l.emplace_back(std::vector<Vector3D>({A2,B2,C2}));
    }
    l.clear();
    l = nw_l;
}

std::vector<std::vector<Vector3D>> clipVistrum(const Vector3D &A ,const Vector3D &B,const Vector3D &C,const Figuur &vistrum) {
    std::vector<std::vector<Vector3D>> l;
    l.emplace_back(std::vector<Vector3D>({A,B,C}));
    for (int i = 0; i <6; i++){
        clipPane(vistrum.vlakInPoints(i), i,l);
    }
    return l;
}
#endif //ENGINE_CLIPPING_H
