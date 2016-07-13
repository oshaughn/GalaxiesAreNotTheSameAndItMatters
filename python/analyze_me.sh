
HERE=`pwd`

# BHBH
mkdir bh_bh; cd bh_bh
for my_type  in exciting boring
do
  mkdir ${my_type}; cd ${my_type}
  ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/${my_type}-galaxy/z0.016.starmetal.dat --type-bbh > my.log
  cd ..
done
mkdir dwarf1; cd dwarf1
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/Dwarf1/dwarf1.z0.starmetal.dat --type-bbh > my.log
   cd ..
mkdir dwarf2; cd dwarf2
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/Dwarf2/dwarf2.z0.1.starmetal.dat --type-bbh > my.log
   cd ..
 mkdir h239; cd h239
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/h239/z0.starmetal.dat --type-bbh > my.log
 cd ..
 mkdir h285; cd h285
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/h285/z0.starmetal.dat --type-bbh > my.log
 cd ..
cd ..

# BHNS
mkdir bh_ns; cd bh_ns
for my_type  in exciting boring
do
  mkdir ${my_type}; cd ${my_type}
  ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/${my_type}-galaxy/z0.016.starmetal.dat --type-bhns > my.log
  cd ..
done
mkdir dwarf1; cd dwarf1
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/Dwarf1/dwarf1.z0.starmetal.dat --type-bhns > my.log
   cd ..
mkdir dwarf2; cd dwarf2
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/Dwarf2/dwarf2.z0.1.starmetal.dat --type-bhns > my.log
   cd ..
 mkdir h239; cd h239
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/h239/z0.starmetal.dat --type-bhns > my.log
 cd ..
 mkdir h285; cd h285
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/h285/z0.starmetal.dat --type-bhns > my.log
 cd ..

cd ..


# BHNS
mkdir ns_ns; cd ns_ns
for my_type  in exciting boring
do
  mkdir ${my_type}; cd ${my_type}
  ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/${my_type}-galaxy/z0.016.starmetal.dat  > my.log
  cd ..
done
mkdir dwarf1; cd dwarf1
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/Dwarf1/dwarf1.z0.starmetal.dat  > my.log
   cd ..
mkdir dwarf2; cd dwarf2
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/Dwarf2/dwarf2.z0.1.starmetal.dat  > my.log
   cd ..
 mkdir h239; cd h239
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/h239/z0.starmetal.dat > my.log
 cd ..
 mkdir h285; cd h285
    ln -sf ${HERE}/mass_max_of_z.dat .
  python ${HERE}/sfr_to_rate_nobins_BBH.py --starmetal-file ~/PersonalJBArchive/h285/z0.starmetal.dat  > my.log
 cd ..
cd ..
