#include <iostream>
//#include <map>
#include <string>
#include <fstream>
#include <unordered_map>
#include <stdio.h>
std::unordered_map<std::string,double> mymap;

extern "C" 
{
  void set_tbl_();
  double lookforkey_(int *, int *, int *, int *, int *,int *, int *, int *, int *, int *);
}

double lookforkey_(int *l1, int *m1, int *l2, int *m2, int *l3, int *m3, int *n1, int *n2, int *n3, int *n4)
{
  std::string sl1 = std::to_string(*l1)+" "+std::to_string(*m1)+" ";
  std::string sl2 = std::to_string(*l2)+" "+std::to_string(*m2)+" ";
  std::string sl3 = std::to_string(*l3)+" "+std::to_string(*m3)+" ";
  std::string sl4 = std::to_string(*n1)+" "+std::to_string(*n2)+" ";
  std::string sl5 = std::to_string(*n3)+" "+std::to_string(*n4);
 
  std::string key = sl1+sl2+sl3+sl4+sl5;
  std::unordered_map<std::string,double>::iterator it;

  it = mymap.find(key);
  if (it == mymap.end()) {
    return 0.0;
  }
  return it->second;
}

void set_tbl_()
{
  std::string line;
  std::ifstream myfile("/afs/mpa-garching.mpg.de/data/haakoan/gw_int/sph_int/table.dat");
  int l1,m1,l2,m2,l3,m3,n1,n2,n3,n4;
  double val;
  if (myfile.is_open())
  {
    while (myfile>>l1>>m1>>l2>>m2>>l3>>m3>>n1>>n2>>n3>>n4>>val )
    {
      std::string sl1 = std::to_string(l1)+" "+std::to_string(m1)+" ";
      std::string sl2 = std::to_string(l2)+" "+std::to_string(m2)+" ";
      std::string sl3 = std::to_string(l3)+" "+std::to_string(m3)+" ";
      std::string sl4 = std::to_string(n1)+" "+std::to_string(n2)+" ";
      std::string sl5 = std::to_string(n3)+" "+std::to_string(n4);
      std::string key = sl1+sl2+sl3+sl4+sl5;
      mymap[key] = val;
      //std::cout << std::fixed << val << "    " << key << '\n';
    }
    myfile.close();
  }
  return;
}

//-lstdc++







  // igr = 0.0d0
  // switch(n3) {
  // case 0:
  // switch(n4) {
  // case 0:
  // igr = 0.0d0;
  // break;
  // case 1:
  //       if((m1+m2+m3 .eq. -1) .or. (m1+m2+m3 .eq. 1)) then
  //          ft = pc_pi
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  // case 2:
  //       if((m1+m2+m3 .eq. -2) .or. (m1+m2+m3 .eq. 2)) then
  //          ft = pc_pi*0.5d0
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //       if((m1+m2+m3 .eq. 0)) then
  //          ft = pc_pi
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //    end select
  // case 1:
  // igr = 0d0;
  // case 2:
  //    switch(n4)
  //    case(0)
  //       if((m1+m2+m3 .eq. -2) .or. (m1+m2+m3 .eq. 2)) then
  //          ft = -pc_pi*0.5d0
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //       if((m1+m2+m3 .eq. 0)) then
  //          ft = pc_pi
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //    case(1)
  //       if((m1+m2+m3 .eq. -3) .or. (m1+m2+m3 .eq. 3)) then
  //          ft = -pc_pi*0.25d0
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //       if((m1+m2+m3 .eq. -1) .or. (m1+m2+m3 .eq. 1)) then
  //          ft = pc_pi*0.25d0
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //    case(2)
  //       if((m1+m2+m3 .eq. -4) .or. (m1+m2+m3 .eq. 4)) then
  //          ft = -pc_pi*0.125d0
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
  //       if((m1+m2+m3 .eq. 0)) then
  //          ft = pc_pi*0.25d0
  //          igr = I3(l1,m1,l2,m2,l3,m3,n1,n2)*ft
  //       end if
 
  // if((*m1+*m2+*m3+*l1+*l2+*l3+*n1)%2 == 0) {return 0.0;}
  // if(*n3 == 0 && *n4 == 0) {return 0.0;}
  // if(*n3 == 0 && *n4 == 2 && ((*m1+*m2+*m3 == -2) || (*m1+*m2+*m3 == 2) || (*m1+*m2+*m3 == 0) )) {
  //   it = mymap.find(key); }
  // if(*n3 == 0 && *n4 == 1 && ((*m1+*m2+*m3 == 1) || (*m1+*m2+*m3 == -1))) {
  //   it = mymap.find(key); }
  // if(*n3 == 1) {return 0.0;}

  // if(*n3 == 2 && *n4 == 0 && ((*m1+*m2+*m3 == -2) || (*m1+*m2+*m3 == 2) || (*m1+*m2+*m3 == 0) )) {
  //   it = mymap.find(key); }
  // if(*n3 == 2 && *n4 == 1 && ((*m1+*m2+*m3 == 1) || (*m1+*m2+*m3 == -1) || (*m1+*m2+*m3 == 3) || (*m1+*m2+*m3 == -3))) {
  //   it = mymap.find(key); }
  // if(*n3 == 2 && *n4 == 2 && ((*m1+*m2+*m3 == -4) || (*m1+*m2+*m3 == 4))) {
  //   it = mymap.find(key); }



  
  // switch(*n3) {
  // case 0:
  //   switch(*n4) {
  //   case 0:
  //     return 0.0;
  //     break;
  //   case 1:
  //     if((*m1+*m2+*m3 == -1) || (*m1+*m2+*m3 == 1)) {
  // 	it = mymap.find(key);}
  //     break;
  //   case 2:
  //     if((*m1+*m2+*m3 == -2) || (*m1+*m2+*m3 == 2) || (*m1+*m2+*m3 == 0) ) {
  // 	it = mymap.find(key);}
  //     break;
  //   }
  //   break;
  // case 1:
  //   return 0.0;
  //   break;
  // case 2:
  //   switch(*n4) {
  //   case 0:
  //     if((*m1+*m2+*m3 == -2) || (*m1+*m2+*m3 == 2) || (*m1+*m2+*m3 == 0) ) {
  // 	it = mymap.find(key);}
  //     break;
  //   case 1:
  //     if((*m1+*m2+*m3 == -3) || (*m1+*m2+*m3 == 3) || (*m1+*m2+*m3 == 1) || (*m1+*m2+*m3 == -1) ) {
  // 	it = mymap.find(key);}
  //     break;
  //   case 2:
  //     if((*m1+*m2+*m3 == -4) || (*m1+*m2+*m3 == 4) || (*m1+*m2+*m3 == 0) ) {
  // 	it = mymap.find(key);}
  //     break;
  //   }
  //   break;
  // }
