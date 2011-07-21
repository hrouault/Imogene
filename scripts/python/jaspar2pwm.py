#!/usr/bin/python

import numpy as N

fmotifs=open(sys.argv[1])

motifs=fmotifs.readlines()

   void
jaspar2pwm()
{

   ifstream fmotifs;
   fmotifs.open(args_info.jaspar2pwm_arg);
   ofstream out;
   out.open("jaspar2pwm-output-wname.dat");
   string dum;
   getline(fmotifs,dum);
   while (!fmotifs.eof()){
      Motif mot;
      mot.name=dum.substr(1);
      mot.id=mot.name.substr(0,mot.name.find(" "));
      while (mot.name.find(" ")!=string::npos){
         mot.name.replace(mot.name.find(" "),1,"_");
      }
      //cout << mot.id << " " << mot.name << endl;

      vd rowdum;
      vvd mat(4,rowdum);
      // read the matrix
      for (int j=0;j<4;j++){
         getline(fmotifs,dum);
         stringstream row(dum);
         row >> dum; // base
         row >> dum; // [
         if (dum.size()>1){ // in case first number is stuck to the [
            dum=dum.substr(1);
         } else {
            row >> dum; // freq
         }
         while (dum.find(']')==string::npos){ // until ]
            stringstream dumtonum(dum);
            double num;
            dumtonum >> num;
            mat[j].push_back(num);
            row >> dum;
         }
         //cout << mat[j] ;
      }

      //cout << "Computing pseudo-count..." << endl;
      mot.motwidth=mat[0].size();
      width=mot.motwidth;
      double sc=scorethr2;
      scorethr2=sc/10*width;
      compalpha();
      scorethr2=sc;

      vd rowzeros(4,0.);
      mot.matprec=vvd(mot.motwidth,rowzeros);
      for (int i=0;i<mot.motwidth;i++){
         double sum=0;
         sum+=mat[0][i];
         sum+=mat[1][i];
         sum+=mat[2][i];
         sum+=mat[3][i];
         sum+=2*alpha+2*beta;
         mot.matprec[i][0]=log((mat[0][i]+alpha)/sum/conca);
         mot.matprec[i][1]=log((mat[3][i]+alpha)/sum/conca);
         mot.matprec[i][2]=log((mat[1][i]+beta)/sum/concc);
         mot.matprec[i][3]=log((mat[2][i]+beta)/sum/concc);

         string letter="A";
         double max(mot.matprec[i][0]);
         if (mot.matprec[i][1]>max) letter="T";
         else if (mot.matprec[i][2]>max) letter="C";
         else if (mot.matprec[i][3]>max) letter="G";
         mot.bsinit.append(letter);
      }
      //cout << mot.matprec << endl;
      mot.displaywname(out);
      getline(fmotifs,dum);
   }
   out.close();
   fmotifs.close();

}

