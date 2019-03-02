module modul1
contains

subroutine naplnitkonstanty(l,xconst) !inicialization of constants depending on spectral parameter l (in thesis denoted j)
real(8), dimension(40), intent(out)  :: xconst
real(8), intent(in) :: l
xconst(1)=sqrt(l) !A
xconst(2)=-sqrt(l)*(l-1)
xconst(3)=-sqrt(l+1)
xconst(4)=-sqrt(l+1)*(l+2)
xconst(5)=-2*sqrt((l-1)/(2*l-1))
xconst(6)=xconst(5)*l !F
xconst(7)=-2*sqrt((l-1)/(2*(2*l+1))) !G
xconst(8)=xconst(7)*(l+1) !H
xconst(9)=2*sqrt((l+1)*(2*l+3)/(6*(2*l-1)*(2*l+1))) !I
xconst(10)=-xconst(9)*(l-1) !J
xconst(11)=-2*sqrt(l*(2*l-1)/(6*(2*l+1)*(2*l+3))) !K
xconst(12)=xconst(11)*(l+2)
xconst(13)=2*sqrt((l+2)/(2*(2*l+1)))
xconst(14)=-xconst(13)*l
xconst(15)=2*sqrt((l+2)/(2*l+3)) !O
xconst(16)=-xconst(15)*(l+1)
xconst(17)=-sqrt(l/(3*(2*l+1)))
xconst(18)=xconst(17)*(l+1)
xconst(19)=sqrt((l-1)/(2*l-1))
xconst(20)=-xconst(19)*(l-2) !T
xconst(21)=-sqrt((l+1)*(2*l+3)/(6*(2*l-1)*(2*l+1)))
xconst(22)=xconst(21)*(l+1)
xconst(23)=-sqrt((l-1)/(2*(2*l+1))) !W
xconst(24)=-xconst(23)*(l-1)
xconst(25)=-sqrt((l+2)/(2*(2*l+1)))
xconst(26)=xconst(25)*(l+2) !Z
xconst(27)=sqrt((l+1)/(3*(2*l+1))) !alpha
xconst(28)=-xconst(27)*l
xconst(29)=sqrt((l*(2*l-1))/(6*(2*l+1)*(2*l+3)))
xconst(30)=-xconst(29)*l
xconst(31)=-sqrt((l+2)/(2*l+3))
xconst(32)=xconst(31)*(l+3)
xconst(33)=-sqrt(l/(3*(2*l+1))) !B1
xconst(34)=sqrt((l-1)/(2*l-1))
xconst(35)=-sqrt((l+1)*(2*l+3)/(6*(2*l+1)*(2*l-1)))
xconst(36)=sqrt((l+1)/(3*(2*l+1)))
xconst(37)=sqrt((l*(2*l-1))/(6*(2*l+1)*(2*l+3)))
xconst(38)=-sqrt((l+2)/(2*l+3))
xconst(39)=sqrt((l-1)/(2*(2*l+1)))
xconst(40)=-sqrt((l+2)/(2*(2*l+1)))

end subroutine naplnitkonstanty

!calculates ur
subroutine vypocitejur(pocetvrstev,l,b,ur)
integer i,pocetvrstev
real(8), dimension(:) :: b,ur
real(8) l
do i=1, pocetvrstev+1
ur(i)=sqrt(l/(2*l+1))*b(6*i-5)-sqrt((l+1)/(2*l+1))*b(6*i-4)
end do
end subroutine

!calculates power
subroutine vypocitejpcelk(pcelk,pocetvrstev,pocetkroku,p,deltat,eta,r,krokuvperiode)
real(8) pcelk,deltat
integer pocetvrstev,i,pocetkroku,tn,krokuvperiode
real(8), dimension(:) :: r,eta
real(8), dimension(:,:) :: p
pcelk=0
do i=0,pocetvrstev-1
do tn=pocetkroku-krokuvperiode, pocetkroku
pcelk=pcelk+p(i+1,tn)*deltat
enddo
pcelk=pcelk-p(i+1,pocetkroku-krokuvperiode)/2-p(i+1,pocetkroku)/2
p(i+1,1)=pcelk/2/eta(i+1) !storing value
pcelk=0
enddo
do i=0,pocetvrstev-1
pcelk=p(i+1,1)*(r(i+1)+r(i+2))/2*(r(i+1)+r(i+2))/2*(r(i+1)-r(i+2))+pcelk
enddo
pcelk=pcelk-p(1,1)*(r(1)+r(2))/2*(r(1)+r(2))/2/2*(r(1)-r(2))-p(pocetvrstev-1,1)*(r(pocetvrstev+1)+r(pocetvrstev))&
/2*(r(pocetvrstev)+r(pocetvrstev+1))/2/2*(r(pocetvrstev)-r(pocetvrstev+1))
pcelk=pcelk/(pocetkroku-krokuvperiode)/deltat
end subroutine

!implementing tidal potential
subroutine zadanislapovehopotencialu(pocetvrstev,l,b,T)
real(8) l,T
integer pocetvrstev
real(8), dimension(:) :: b
do i=1, pocetvrstev*6+2
b(i)=0
end do
b(pocetvrstev*6+1)=-T*sqrt(l/(2*l+1))
b(pocetvrstev*6+2)=T*sqrt((l+1)/(2*l+1))
end subroutine

!implementing force on the top
subroutine zadanisily(pocetvrstev,l,b,T)
real(8) l,T
integer pocetvrstev
real(8), dimension(:) :: b
do i=1, pocetvrstev*6+2
b(i)=0
end do
b(1)=-T*sqrt(l/(2*l+1))
b(2)=T*sqrt((l+1)/(2*l+1))
end subroutine
end module

program naplnenimatice
use modul1
implicit none
!inicialization
real(8), dimension(:,:), allocatable :: p,eps,epsc,epsnc,epsmemory
integer :: pocetvrstev,i,j,tn,k,n,m,pocetkroku,krokuvperiode,planet,viscositystart,viscosityend
real(8), dimension(:), allocatable :: r,rhat,b,b2,b2nc,b2c,b3,b3c,b3nc,bmemory,bmemory3, rtilde,rhatd,rtilded,ur,mu,eta,ur22
real(8), dimension(:), allocatable :: bcmemory, bncmemory
real(8) :: l,rmax,rmin,rpomoc,g,g2,T,T2,T2nc,T2c,e,rho,deltarho,rhocore,omega,mukonst,deltat,pi
real(8), dimension(40) :: xconst
real(8) :: d,pcelk, urmax, tnmax
integer, dimension(:), allocatable :: indx
real(8), dimension(:,:), allocatable :: al,a,a2,a3
real(8) :: urmemory, urmaxmemory, theta

l=2 !spectral parameter l (in thesis denoted j)
m=0 !spectral parameter m
pi=3.1415926535897932384626433832795d0 !pi constant
krokuvperiode=2000.0d0 !numerical discretization in time, number of steps in one period
pocetkroku=10*krokuvperiode !number of total steps in time

print*, "Choose from examined planetary bodies"
print*, "Enceladus 26km(1), Enceladus 52km(2), Europa(3), Exoplanet 1day(4)"
print*, "Exoplaneta 5days(5), Exoplaneta 20days(6), Exoplaneta 50days(7), Merkur(8)"
read*, planet

SELECT CASE (planet)
CASE (1:2)!Enceladus
g=0.11d0!gravitational acceleration at surface
g2=0.13d0!gravitational acceleration at the interface between mantle and core
rho=925d0!density of mantle
deltarho=82d0!difference between the core density and mantle density
rhocore=rho+deltarho
T=1000*rho*g!testing force
e=0.0045d0!orbital eccentricity
omega=2*pi/(1.370218d0*60*60*24)!orbital angular frequence
rmax=252100!top boundary in meters
pocetvrstev=20!number of layers in vertical spatial discretization
if (planet==2) then
rmin=rmax-52000!52km
else
rmin=rmax-26000!26km
endif
mukonst=3.3d9!shear modulus
viscositystart=10!beginning of viscosity measurement
viscosityend=22!end of viscosity measurement
CASE (3)!Europa
g=1.314d0
g2=1.314d0
rho=925d0
deltarho=82d0
rhocore=rho+deltarho
T=1000*rho*g
e=0.009d0
omega=2*pi/(3.551d0*60*60*24)
rmax=1561000
pocetvrstev=20
rmin=rmax-30000
mukonst=3.3d9
viscositystart=10
viscosityend=22
CASE (4:7)!Exoplanet
g=9.72868d0
g2=10.73510d0
rho=4500d0
deltarho=7500d0
rhocore=rho+deltarho
T=1000*rho*g
e=0.1d0

rmax=6400000
pocetvrstev=20
rmin=rmax-3200000
mukonst=70.0d9
viscositystart=14
viscosityend=24
SELECT CASE (planet)!different orbital periods
CASE(4)
omega=2*pi/(1.0d0*60*60*24)
CASE(5)
omega=2*pi/(5.0d0*60*60*24)
CASE(6)
omega=2*pi/(20.0d0*60*60*24)
CASE(7)
omega=2*pi/(50.0d0*60*60*24)
END SELECT

CASE (8)!mercury
g=3.7d0
g2=3.7d0
rho=4000d0
deltarho=3000d0
rhocore=rho+deltarho
T=1000*rho*g
e=0.2d0
omega=2*pi/(87.969d0*60*60*24)
rmax=2450000
pocetvrstev=20
rmin=rmax-400000
mukonst=70.0d9
viscositystart=14
viscosityend=24
END SELECT


call naplnitkonstanty(l,xconst)

!Elasticity part
allocate(r(pocetvrstev+1))
allocate(rhat(pocetvrstev))
allocate(rtilde(pocetvrstev))
allocate(mu(pocetvrstev))
allocate(rhatd(pocetvrstev-1))
allocate(rtilded(pocetvrstev-1))
allocate(eta(pocetvrstev))

!define auxiliary variables in order to simplify matrix A as done in thesis
do i=1, pocetvrstev+1
r(i)=rmax-(rmax-rmin)/((pocetvrstev-1)*2)*(2*i-3)
enddo
do i=1, pocetvrstev
rhat(i)=r(i+1)+r(i)
rtilde(i)=r(i+1)-r(i)
enddo
do i=1, pocetvrstev-1
rhatd(i)=r(i+2)+r(i)
rtilded(i)=(r(i+2)-r(i))/2
enddo

!fill shear modulus constant in all layers
do i=1, pocetvrstev
mu(i)=mukonst
end do

allocate(indx(pocetvrstev*6+2))
allocate(a(pocetvrstev*6+2,15))
allocate(al(pocetvrstev*6+2,7))
allocate(b(pocetvrstev*6+2))
allocate(ur(pocetvrstev+1))
allocate(ur22(pocetvrstev+1))

!right hand side external force on the top boundary, possibility of imposing force by changing last parameter
call zadanisily(pocetvrstev,l,b,0d0)
a=0!inicialization of matrix in bandform

do i=0, pocetvrstev-1
!continuity equation
a(6*i+2+1,4-1+3)=-rhat(i+1)*xconst(1)+rtilde(i+1)*xconst(2)
a(6*i+2+1,4+0+3)=-rhat(i+1)*xconst(3)+rtilde(i+1)*xconst(4)
a(6*i+2+1,4+5+3)=rhat(i+1)*xconst(1)+rtilde(i+1)*xconst(2)
a(6*i+2+1,4+6+3)=rhat(i+1)*xconst(3)+rtilde(i+1)*xconst(4)
!rheological relationship
a(6*i+2+2,3-1+3)=(-rhat(i+1)*mu(i+1)*xconst(5)+rtilde(i+1)*mu(i+1)*xconst(6))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+2,3+5+3)=(rhat(i+1)*mu(i+1)*xconst(5)+rtilde(i+1)*mu(i+1)*xconst(6))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+2,3+1+3)=(rhat(i+1)*rtilde(i+1))/(rhat(i+1)*rtilde(i+1))

a(6*i+2+3,2-1+3)=(-rhat(i+1)*mu(i+1)*xconst(9)+rtilde(i+1)*mu(i+1)*xconst(10))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+3,2+0+3)=(-rhat(i+1)*mu(i+1)*xconst(11)+rtilde(i+1)*mu(i+1)*xconst(12))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+3,2+5+3)=(rhat(i+1)*mu(i+1)*xconst(9)+rtilde(i+1)*mu(i+1)*xconst(10))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+3,2+6+3)=(rhat(i+1)*mu(i+1)*xconst(11)+rtilde(i+1)*mu(i+1)*xconst(12))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+3,2+2+3)=(rhat(i+1)*rtilde(i+1))/(rhat(i+1)*rtilde(i+1))

a(6*i+2+4,1+0+3)=(-rhat(i+1)*mu(i+1)*xconst(15)+rtilde(i+1)*mu(i+1)*xconst(16))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+4,1+6+3)=(rhat(i+1)*mu(i+1)*xconst(15)+rtilde(i+1)*mu(i+1)*xconst(16))/(rhat(i+1)*rtilde(i+1))
a(6*i+2+4,1+4+3)=(rhat(i+1)*rtilde(i+1))/(rhat(i+1)*rtilde(i+1))
!equationofmotion
if (i<pocetvrstev-1) then
a(6*i+2+5,1+3)=-rhatd(i+1)*xconst(19)+rtilded(i+1)*xconst(20)
a(6*i+2+5,2+3)=-rhatd(i+1)*xconst(21)+rtilded(i+1)*xconst(22)
a(6*i+2+5,3+3)=-rhatd(i+1)*xconst(17)+rtilded(i+1)*xconst(18)

a(6*i+2+5,1+6+3)=rhatd(i+1)*xconst(19)+rtilded(i+1)*xconst(20)
a(6*i+2+5,2+6+3)=rhatd(i+1)*xconst(21)+rtilded(i+1)*xconst(22)
a(6*i+2+5,3+6+3)=rhatd(i+1)*xconst(17)+rtilded(i+1)*xconst(18)
end if

if (i<pocetvrstev-1) then
a(6*i+2+6,-1+2+3)=-rhatd(i+1)*xconst(29)+rtilded(i+1)*xconst(30)
a(6*i+2+6,-1+3+3)=-rhatd(i+1)*xconst(27)+rtilded(i+1)*xconst(28)
a(6*i+2+6,-1+4+3)=-rhatd(i+1)*xconst(31)+rtilded(i+1)*xconst(32)

a(6*i+2+6,-1+2+6+3)=rhatd(i+1)*xconst(29)+rtilded(i+1)*xconst(30)
a(6*i+2+6,-1+3+6+3)=rhatd(i+1)*xconst(27)+rtilded(i+1)*xconst(28)
a(6*i+2+6,-1+4+6+3)=rhatd(i+1)*xconst(31)+rtilded(i+1)*xconst(32)
end if
end do

!top boundary condition
a(1,5+3)=rho*g*l/((2*l+1)*2)
a(1,6+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a(1,7+3)=xconst(34)
a(1,8+3)=xconst(35)
a(1,9+3)=xconst(33)
a(1,11+3)=rho*g*l/((2*l+1)*2)
a(1,12+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))

a(2,7+3)=xconst(37)
a(2,8+3)=xconst(36)
a(2,9+3)=xconst(38)
a(2,4+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a(2,5+3)=rho*g*(l+1)/((2*l+1)*2)
a(2,10+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a(2,11+3)=rho*g*(l+1)/((2*l+1)*2)

!bottom boundary condition
a(pocetvrstev*6+1,1+3)=xconst(34)
a(pocetvrstev*6+1,2+3)=xconst(35)
a(pocetvrstev*6+1,3+3)=xconst(33)
a(pocetvrstev*6+1,-1+3)=-deltarho*g2*l/((2*l+1)*2)
a(pocetvrstev*6+1,0+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a(pocetvrstev*6+1,5+3)=-deltarho*g2*l/((2*l+1)*2)
a(pocetvrstev*6+1,6+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))

a(pocetvrstev*6+2,1+3)=xconst(37)
a(pocetvrstev*6+2,2+3)=xconst(36)
a(pocetvrstev*6+2,3+3)=xconst(38)
a(pocetvrstev*6+2,-2+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a(pocetvrstev*6+2,-1+3)=-deltarho*g2*(l+1)/((2*l+1)*2)
a(pocetvrstev*6+2,4+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a(pocetvrstev*6+2,5+3)=-deltarho*g2*(l+1)/((2*l+1)*2)

call bandecdp(a,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,d)!procedure from numerical recipes to solve matrix in band form
call banbksdp(a,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b)!procedure from numerical recipes to calculate matrix equation Ax=b

!----------------------------------------------------------------------------------------------------------
!End of model of elasticity
!----------------------------------------------------------------------------------------------------------

allocate(b2(pocetvrstev*6+2))
allocate(bmemory(pocetvrstev*6+2))
allocate(a2(pocetvrstev*6+2,15))
allocate(p(pocetvrstev,pocetkroku))
allocate(b2c(pocetvrstev*6+2))
allocate(b2nc(pocetvrstev*6+2))
allocate(bcmemory(pocetvrstev*6+2))
allocate(bncmemory(pocetvrstev*6+2))
allocate(b3(pocetvrstev*6+2))
allocate(bmemory3(pocetvrstev*6+2))
allocate(a3(pocetvrstev*6+2,15))
allocate(eps(pocetvrstev*6+2,pocetkroku+1))
allocate(b3c(pocetvrstev*6+2))
allocate(b3nc(pocetvrstev*6+2))
allocate(epsmemory(pocetvrstev*6+2,pocetkroku+1))
allocate(epsc(pocetvrstev*6+2,pocetkroku+1))
allocate(epsnc(pocetvrstev*6+2,pocetkroku+1))

open(1, file="pcelketa")!file where power and viscosity is stored for maxwell model
open(2, file="pcelketakelvin") !file where power and viscosity is stored for kelvin model
open(4, file="ureta")!file where radial displacement(ur_20) and viscosity is stored for maxwell model
open(5, file="uretakelvin")!file where radial displacement(ur_20) and viscosity is stored for kelvin model
open(11, file="tneta")!file where phase offset and viscosity is stored for maxwell model
open(12, file="tnetakelvin")!file where phase offset and viscosity is stored for kelvin model
open(17, file="urtheta")!file where power and viscosity is stored for maxwell model
open(18, file="urthetakelvin")!file where power and viscosity is stored for kelvin model

do j=viscositystart*20,viscosityend*20  !viscosity 10-22 ice, 14-24 silicate
    eta=10.0d0**(j/20.0d0)

deltat=2*pi/omega/krokuvperiode!timestepchosensmall

a2=0
do i=0, pocetvrstev-1
!continuity equation
a2(6*i+2+1,4-1+3)=-rhat(i+1)*xconst(1)+rtilde(i+1)*xconst(2)
a2(6*i+2+1,4+0+3)=-rhat(i+1)*xconst(3)+rtilde(i+1)*xconst(4)
a2(6*i+2+1,4+5+3)=rhat(i+1)*xconst(1)+rtilde(i+1)*xconst(2)
a2(6*i+2+1,4+6+3)=rhat(i+1)*xconst(3)+rtilde(i+1)*xconst(4)
!rheological relationship
a2(6*i+2+2,3-1+3)=(-rhat(i+1)*mu(i+1)*xconst(5)+rtilde(i+1)*mu(i+1)*xconst(6))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+2,3+5+3)=(rhat(i+1)*mu(i+1)*xconst(5)+rtilde(i+1)*mu(i+1)*xconst(6))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+2,3+1+3)=(rhat(i+1)*rtilde(i+1)*(1+mu(i+1)/eta(i+1)*deltat/2.0D0))/(rhat(i+1)*rtilde(i+1))

a2(6*i+2+3,2-1+3)=(-rhat(i+1)*mu(i+1)*xconst(9)+rtilde(i+1)*mu(i+1)*xconst(10))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+3,2+0+3)=(-rhat(i+1)*mu(i+1)*xconst(11)+rtilde(i+1)*mu(i+1)*xconst(12))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+3,2+5+3)=(rhat(i+1)*mu(i+1)*xconst(9)+rtilde(i+1)*mu(i+1)*xconst(10))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+3,2+6+3)=(rhat(i+1)*mu(i+1)*xconst(11)+rtilde(i+1)*mu(i+1)*xconst(12))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+3,2+2+3)=(rhat(i+1)*rtilde(i+1)*(1+mu(i+1)/eta(i+1)*deltat/2.0D0))/(rhat(i+1)*rtilde(i+1))

a2(6*i+2+4,1+0+3)=(-rhat(i+1)*mu(i+1)*xconst(15)+rtilde(i+1)*mu(i+1)*xconst(16))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+4,1+6+3)=(rhat(i+1)*mu(i+1)*xconst(15)+rtilde(i+1)*mu(i+1)*xconst(16))/(rhat(i+1)*rtilde(i+1))
a2(6*i+2+4,1+4+3)=(rhat(i+1)*rtilde(i+1)*(1+mu(i+1)/eta(i+1)*deltat/2.0D0))/(rhat(i+1)*rtilde(i+1))


!equation of motion
if (i<pocetvrstev-1) then
a2(6*i+2+5,1+3)=(-rhatd(i+1)*xconst(19)+rtilded(i+1)*xconst(20))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+5,2+3)=(-rhatd(i+1)*xconst(21)+rtilded(i+1)*xconst(22))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+5,3+3)=(-rhatd(i+1)*xconst(17)+rtilded(i+1)*xconst(18))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+5,1+6+3)=(rhatd(i+1)*xconst(19)+rtilded(i+1)*xconst(20))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+5,2+6+3)=(rhatd(i+1)*xconst(21)+rtilded(i+1)*xconst(22))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+5,3+6+3)=(rhatd(i+1)*xconst(17)+rtilded(i+1)*xconst(18))/(rhatd(i+1)*rtilded(i+1))

end if

if (i<pocetvrstev-1) then
a2(6*i+2+6,-1+2+3)=(-rhatd(i+1)*xconst(29)+rtilded(i+1)*xconst(30))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+6,-1+3+3)=(-rhatd(i+1)*xconst(27)+rtilded(i+1)*xconst(28))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+6,-1+4+3)=(-rhatd(i+1)*xconst(31)+rtilded(i+1)*xconst(32))/(rhatd(i+1)*rtilded(i+1))

a2(6*i+2+6,-1+2+6+3)=(rhatd(i+1)*xconst(29)+rtilded(i+1)*xconst(30))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+6,-1+3+6+3)=(rhatd(i+1)*xconst(27)+rtilded(i+1)*xconst(28))/(rhatd(i+1)*rtilded(i+1))
a2(6*i+2+6,-1+4+6+3)=(rhatd(i+1)*xconst(31)+rtilded(i+1)*xconst(32))/(rhatd(i+1)*rtilded(i+1))

end if
end do

!boundary condition
a2(1,5+3)=rho*g*l/((2*l+1)*2)
a2(1,6+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a2(1,7+3)=xconst(34)
a2(1,8+3)=xconst(35)
a2(1,9+3)=xconst(33)
a2(1,11+3)=rho*g*l/((2*l+1)*2)
a2(1,12+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))


a2(2,7+3)=xconst(37)
a2(2,8+3)=xconst(36)
a2(2,9+3)=xconst(38)
a2(2,4+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a2(2,5+3)=rho*g*(l+1)/((2*l+1)*2)
a2(2,10+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a2(2,11+3)=rho*g*(l+1)/((2*l+1)*2)

!bottom boundary condition
a2(pocetvrstev*6+1,1+3)=xconst(34)
a2(pocetvrstev*6+1,2+3)=xconst(35)
a2(pocetvrstev*6+1,3+3)=xconst(33)
a2(pocetvrstev*6+1,-1+3)=-deltarho*g2*l/((2*l+1)*2)
a2(pocetvrstev*6+1,0+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a2(pocetvrstev*6+1,5+3)=-deltarho*g2*l/((2*l+1)*2)
a2(pocetvrstev*6+1,6+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))


a2(pocetvrstev*6+2,1+3)=xconst(37)
a2(pocetvrstev*6+2,2+3)=xconst(36)
a2(pocetvrstev*6+2,3+3)=xconst(38)
a2(pocetvrstev*6+2,-2+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a2(pocetvrstev*6+2,-1+3)=-deltarho*g2*(l+1)/((2*l+1)*2)
a2(pocetvrstev*6+2,4+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a2(pocetvrstev*6+2,5+3)=-deltarho*g2*(l+1)/((2*l+1)*2)

call bandecdp(a2,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,d)

p=0
bmemory=0
bncmemory=0
bcmemory=0
urmax=0
urmaxmemory=0

open(3, file='time')!file where radial displacement and time is stored
do tn=1,pocetkroku
!tidal potential for mercury
if (planet.eq.8) then
T2=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2.0d0*(r(pocetvrstev+1)+&
r(pocetvrstev))/2.0d0*omega*omega*e*(-sqrt(9*pi/5.0d0)*(cos(omega*tn*deltat)+3.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
T2nc=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2.0d0*(r(pocetvrstev+1)+&
r(pocetvrstev))/2.0d0*omega*omega*(sqrt(3*pi/10.0d0)*((1+6*e*e)*cos(omega*tn*deltat)-1.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
T2c=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2.0d0*(r(pocetvrstev+1)+&
r(pocetvrstev))/2.0d0*omega*omega*(sqrt(3*pi/10.0d0)*((1-11*e*e)*sin(omega*tn*deltat)-1.0d0/2.0d0*e*sin(2*omega*tn*deltat)))

call zadanislapovehopotencialu(pocetvrstev,l,b2,T2)
call zadanislapovehopotencialu(pocetvrstev,l,b2nc,T2nc)
call zadanislapovehopotencialu(pocetvrstev,l,b2c,T2c)
!!!
do i=1, pocetvrstev-1
b2(1+6*i)=-rho*r(i+1)*omega*omega*e*(-sqrt(18*pi)*(cos(omega*tn*deltat)+3.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
b2nc(1+6*i)=-rho*r(i+1)*omega*omega*(sqrt(3*pi)*((1+6*e*e)*cos(omega*tn*deltat)-1.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
b2c(1+6*i)=-rho*r(i+1)*omega*omega*(sqrt(3*pi)*((1-11*e*e)*sin(omega*tn*deltat)-1.0d0/2.0d0*e*sin(2*omega*tn*deltat)))
end do
else
!tidal potential force in 20, 22, 2-2 decomposed into real and compex part
T2=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2*(r(pocetvrstev+1)+&
r(pocetvrstev))/2*omega*omega*e*(-sqrt(9*pi/5)*cos(omega*tn*deltat))
T2nc=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2*(r(pocetvrstev+1)+&
r(pocetvrstev))/2*omega*omega*e*(sqrt(27*pi/10)*cos(omega*tn*deltat))
T2c=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2*(r(pocetvrstev+1)+&
r(pocetvrstev))/2*omega*omega*e*(-sqrt(24*pi/5)*sin(omega*tn*deltat))

call zadanislapovehopotencialu(pocetvrstev,l,b2,T2)
call zadanislapovehopotencialu(pocetvrstev,l,b2nc,T2nc)
call zadanislapovehopotencialu(pocetvrstev,l,b2c,T2c)

!adding tidal potential to equation of motion
do i=1, pocetvrstev-1
b2(1+6*i)=-rho*r(i+1)*omega*omega*e*(-sqrt(18*pi)*cos(omega*tn*deltat))
b2nc(1+6*i)=-rho*r(i+1)*omega*omega*e*(sqrt(27*pi)*cos(omega*tn*deltat))
b2c(1+6*i)=-rho*r(i+1)*omega*omega*e*(-sqrt(48*pi)*sin(omega*tn*deltat))
end do
endif

!computing maxwell right hand side
do i=0,pocetvrstev-1
b2(4+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(3+6*i)/2+bmemory(3+6*i))
b2(5+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(4+6*i)/2+bmemory(4+6*i))
b2(6+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(6+6*i)/2+bmemory(6+6*i))
b2nc(4+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(3+6*i)/2+bncmemory(3+6*i))
b2nc(5+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(4+6*i)/2+bncmemory(4+6*i))
b2nc(6+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(6+6*i)/2+bncmemory(6+6*i))
b2c(4+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(3+6*i)/2+bcmemory(3+6*i))
b2c(5+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(4+6*i)/2+bcmemory(4+6*i))
b2c(6+6*i)=-mu(i+1)/eta(i+1)*deltat*(b(6+6*i)/2+bcmemory(6+6*i))
enddo

!solving equation
call banbksdp(a2,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b2)
call banbksdp(a2,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b2nc)
call banbksdp(a2,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b2c)

!storing information about the previous timestep
do i=0, pocetvrstev-1
bmemory(3+6*i)=bmemory(3+6*i)+b2(3+6*i)
bmemory(4+6*i)=bmemory(4+6*i)+b2(4+6*i)
bmemory(6+6*i)=bmemory(6+6*i)+b2(6+6*i)

bncmemory(3+6*i)=bncmemory(3+6*i)+b2nc(3+6*i)
bncmemory(4+6*i)=bncmemory(4+6*i)+b2nc(4+6*i)
bncmemory(6+6*i)=bncmemory(6+6*i)+b2nc(6+6*i)

bcmemory(3+6*i)=bcmemory(3+6*i)+b2c(3+6*i)
bcmemory(4+6*i)=bcmemory(4+6*i)+b2c(4+6*i)
bcmemory(6+6*i)=bcmemory(6+6*i)+b2c(6+6*i)

!calculating power
p(i+1,tn)=(b2(3+6*i)*b2(3+6*i))+(b2(4+6*i)*b2(4+6*i))+(b2(6+6*i)*b2(6+6*i))+&
2*((b2nc(3+6*i)*b2nc(3+6*i))+(b2nc(4+6*i)*b2nc(4+6*i))+(b2nc(6+6*i)*b2nc(6+6*i)))+&
2*(b2c(3+6*i)*b2c(3+6*i)+b2c(4+6*i)*b2c(4+6*i)+b2c(6+6*i)*b2c(6+6*i))+&
p(i+1,tn)
end do

call vypocitejur(pocetvrstev,l,b2,ur)
call vypocitejur(pocetvrstev,l,b2nc,ur22)

!calculating precise radial displacement
theta = pi/2
urmemory = (ur(1)+ur(2))/2*1/4*sqrt(5.0d0/pi)*(3*cos(theta)*cos(theta)-1)+&
(ur22(1)+ur22(2))*1/4*sqrt(15.0d0/2/pi)*sin(theta)*sin(theta)

if (tn>pocetkroku-krokuvperiode) then
if (abs(urmemory) > urmaxmemory) then
urmaxmemory = abs(urmemory)
endif
if ((ur(1)+ur(2))/2<urmax) then
    urmax=(ur(1)+ur(2))/2
    tnmax=1.0d0*(tn-(pocetkroku-krokuvperiode))/krokuvperiode
endif
endif
if (tnmax.eq.1) then
tnmax = 0
endif
write(3,*),tn,(ur(1)+ur(2))/2 !ur_20 in time

enddo !tn cycle
close(3)
write(4,*),eta(1),urmax
write(11,*),eta(1),tnmax
write(17,*),eta(1),urmaxmemory

print*,"viscosity ", eta(1)

!calculating power
call vypocitejpcelk(pcelk,pocetvrstev,pocetkroku,p,deltat,eta,r,krokuvperiode)
print*, pcelk
write(1,*), eta(1),pcelk

!End of Maxwell model
!----------------------------------------------------------------------------------------------------------
!Beginning of Kelvin-Voigt model
!----------------------------------------------------------------------------------------------------------

a3=0
do i=0, pocetvrstev-1
!continuity equation
a3(6*i+2+1,4-1+3)=-rhat(i+1)*xconst(1)+rtilde(i+1)*xconst(2)
a3(6*i+2+1,4+0+3)=-rhat(i+1)*xconst(3)+rtilde(i+1)*xconst(4)
a3(6*i+2+1,4+5+3)=rhat(i+1)*xconst(1)+rtilde(i+1)*xconst(2)
a3(6*i+2+1,4+6+3)=rhat(i+1)*xconst(3)+rtilde(i+1)*xconst(4)
!rheological relationship
a3(6*i+2+2,3-1+3)=(-rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(5)+rtilde(i+1)*&
(mu(i+1)+eta(i+1)/deltat)*xconst(6))/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+2,3+5+3)=(rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(5)+rtilde(i+1)*&
(mu(i+1)+eta(i+1)/deltat)*xconst(6))/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+2,3+1+3)=1

a3(6*i+2+3,2-1+3)=(-rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(9)+rtilde(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(10))&
/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+3,2+0+3)=(-rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(11)+rtilde(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(12))&
/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+3,2+5+3)=(rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(9)+rtilde(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(10))&
/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+3,2+6+3)=(rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(11)+rtilde(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(12))&
/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+3,2+2+3)=1

a3(6*i+2+4,1+0+3)=(-rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(15)+rtilde(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(16))&
/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+4,1+6+3)=(rhat(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(15)+rtilde(i+1)*(mu(i+1)+eta(i+1)/deltat)*xconst(16))&
/(rhat(i+1)*rtilde(i+1))
a3(6*i+2+4,1+4+3)=1


!equation of motion
if (i<pocetvrstev-1) then
a3(6*i+2+5,1+3)=(-rhatd(i+1)*xconst(19)+rtilded(i+1)*xconst(20))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+5,2+3)=(-rhatd(i+1)*xconst(21)+rtilded(i+1)*xconst(22))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+5,3+3)=(-rhatd(i+1)*xconst(17)+rtilded(i+1)*xconst(18))/(rhatd(i+1)*rtilded(i+1))

a3(6*i+2+5,1+6+3)=(rhatd(i+1)*xconst(19)+rtilded(i+1)*xconst(20))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+5,2+6+3)=(rhatd(i+1)*xconst(21)+rtilded(i+1)*xconst(22))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+5,3+6+3)=(rhatd(i+1)*xconst(17)+rtilded(i+1)*xconst(18))/(rhatd(i+1)*rtilded(i+1))

end if

if (i<pocetvrstev-1) then
a3(6*i+2+6,-1+2+3)=(-rhatd(i+1)*xconst(29)+rtilded(i+1)*xconst(30))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+6,-1+3+3)=(-rhatd(i+1)*xconst(27)+rtilded(i+1)*xconst(28))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+6,-1+4+3)=(-rhatd(i+1)*xconst(31)+rtilded(i+1)*xconst(32))/(rhatd(i+1)*rtilded(i+1))

a3(6*i+2+6,-1+2+6+3)=(rhatd(i+1)*xconst(29)+rtilded(i+1)*xconst(30))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+6,-1+3+6+3)=(rhatd(i+1)*xconst(27)+rtilded(i+1)*xconst(28))/(rhatd(i+1)*rtilded(i+1))
a3(6*i+2+6,-1+4+6+3)=(rhatd(i+1)*xconst(31)+rtilded(i+1)*xconst(32))/(rhatd(i+1)*rtilded(i+1))

end if
end do
!boundary conditions
a3(1,5+3)=rho*g*l/((2*l+1)*2)
a3(1,6+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a3(1,7+3)=xconst(34)
a3(1,8+3)=xconst(35)
a3(1,9+3)=xconst(33)
a3(1,11+3)=rho*g*l/((2*l+1)*2)
a3(1,12+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))

a3(2,7+3)=xconst(37)
a3(2,8+3)=xconst(36)
a3(2,9+3)=xconst(38)
a3(2,4+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a3(2,5+3)=rho*g*(l+1)/((2*l+1)*2)
a3(2,10+3)=rho*g*(-sqrt(l*(l+1))/((2*l+1)*2))
a3(2,11+3)=rho*g*(l+1)/((2*l+1)*2)

a3(pocetvrstev*6+1,1+3)=xconst(34)
a3(pocetvrstev*6+1,2+3)=xconst(35)
a3(pocetvrstev*6+1,3+3)=xconst(33)
a3(pocetvrstev*6+1,-1+3)=-deltarho*g2*l/((2*l+1)*2)
a3(pocetvrstev*6+1,0+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a3(pocetvrstev*6+1,5+3)=-deltarho*g2*l/((2*l+1)*2)
a3(pocetvrstev*6+1,6+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))

a3(pocetvrstev*6+2,1+3)=xconst(37)
a3(pocetvrstev*6+2,2+3)=xconst(36)
a3(pocetvrstev*6+2,3+3)=xconst(38)
a3(pocetvrstev*6+2,-2+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a3(pocetvrstev*6+2,-1+3)=-deltarho*g2*(l+1)/((2*l+1)*2)
a3(pocetvrstev*6+2,4+3)=-deltarho*g2*(-sqrt(l*(l+1))/((2*l+1)*2))
a3(pocetvrstev*6+2,5+3)=-deltarho*g2*(l+1)/((2*l+1)*2)


call bandecdp(a3,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,d)

eps=0
b3=0
b3c=0
b3nc=0
epsmemory=0
epsc=0
epsnc=0
p=0
urmax=0
urmaxmemory=0

open(3,file='timekelvin')
do tn=1,pocetkroku

!tidal potential for mercury
if (planet.eq.8) then
T2=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2.0d0*(r(pocetvrstev+1)+&
r(pocetvrstev))/2.0d0*omega*omega*e*(-sqrt(9*pi/5.0d0)*(cos(omega*tn*deltat)+3.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
T2nc=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2.0d0*(r(pocetvrstev+1)+&
r(pocetvrstev))/2.0d0*omega*omega*(sqrt(3*pi/10.0d0)*((1+6*e*e)*cos(omega*tn*deltat)-1.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
T2c=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2.0d0*(r(pocetvrstev+1)+&
r(pocetvrstev))/2.0d0*omega*omega*(sqrt(3*pi/10.0d0)*((1-11*e*e)*sin(omega*tn*deltat)-1.0d0/2.0d0*e*sin(2*omega*tn*deltat)))

call zadanislapovehopotencialu(pocetvrstev,l,b3,T2)
call zadanislapovehopotencialu(pocetvrstev,l,b3nc,T2nc)
call zadanislapovehopotencialu(pocetvrstev,l,b3c,T2c)

do i=1, pocetvrstev-1
b2(1+6*i)=-rho*r(i+1)*omega*omega*e*(-sqrt(18*pi)*(cos(omega*tn*deltat)+3.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
b2nc(1+6*i)=-rho*r(i+1)*omega*omega*(sqrt(3*pi)*((1+6*e*e)*cos(omega*tn*deltat)-1.0d0/2.0d0*e*cos(2*omega*tn*deltat)))
b2c(1+6*i)=-rho*r(i+1)*omega*omega*(sqrt(3*pi)*((1-11*e*e)*sin(omega*tn*deltat)-1.0d0/2.0d0*e*sin(2*omega*tn*deltat)))
end do






else
!tidal potential else
T2=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2*(r(pocetvrstev+1)+&
r(pocetvrstev))/2*omega*omega*e*(-sqrt(9*pi/5)*cos(omega*tn*deltat))
T2nc=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2*(r(pocetvrstev+1)+&
r(pocetvrstev))/2*omega*omega*e*(sqrt(27*pi/10)*cos(omega*tn*deltat))
T2c=rhocore*(r(pocetvrstev+1)+r(pocetvrstev))/2*(r(pocetvrstev+1)+&
r(pocetvrstev))/2*omega*omega*e*(-sqrt(24*pi/5)*sin(omega*tn*deltat))

call zadanislapovehopotencialu(pocetvrstev,l,b3,T2)
call zadanislapovehopotencialu(pocetvrstev,l,b3nc,T2nc)
call zadanislapovehopotencialu(pocetvrstev,l,b3c,T2c)
!adding tidal potential to equation of motion
do i=1, pocetvrstev-1
b3(1+6*i)=-rho*r(i+1)*omega*omega*e*(-sqrt(18*pi)*cos(omega*tn*deltat))
b3nc(1+6*i)=-rho*r(i+1)*omega*omega*e*(sqrt(27*pi)*cos(omega*tn*deltat))
b3c(1+6*i)=-rho*r(i+1)*omega*omega*e*(-sqrt(48*pi)*sin(omega*tn*deltat))
end do
endif

!right hand side kelvin
do i=0,pocetvrstev-1
b3(4+6*i)=-2*eta(i+1)/deltat*eps(3+6*i,tn)
b3(5+6*i)=-2*eta(i+1)/deltat*eps(4+6*i,tn)
b3(6+6*i)=-2*eta(i+1)/deltat*eps(6+6*i,tn)
b3c(4+6*i)=-2*eta(i+1)/deltat*epsc(3+6*i,tn)
b3c(5+6*i)=-2*eta(i+1)/deltat*epsc(4+6*i,tn)
b3c(6+6*i)=-2*eta(i+1)/deltat*epsc(6+6*i,tn)
b3nc(4+6*i)=-2*eta(i+1)/deltat*epsnc(3+6*i,tn)
b3nc(5+6*i)=-2*eta(i+1)/deltat*epsnc(4+6*i,tn)
b3nc(6+6*i)=-2*eta(i+1)/deltat*epsnc(6+6*i,tn)
end do
call banbksdp(a3,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b3)
call banbksdp(a3,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b3nc)
call banbksdp(a3,pocetvrstev*6+2,7,7,pocetvrstev*6+2,15,al,7,indx,b3c)

do i=0, pocetvrstev-1
eps(3+6*i,tn+1)=(eps(3+6*i,tn)*2*eta(i+1)/deltat+b3(3+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
eps(4+6*i,tn+1)=(eps(4+6*i,tn)*2*eta(i+1)/deltat+b3(4+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
eps(6+6*i,tn+1)=(eps(6+6*i,tn)*2*eta(i+1)/deltat+b3(6+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))

epsc(3+6*i,tn+1)=(epsc(3+6*i,tn)*2*eta(i+1)/deltat+b3c(3+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
epsc(4+6*i,tn+1)=(epsc(4+6*i,tn)*2*eta(i+1)/deltat+b3c(4+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
epsc(6+6*i,tn+1)=(epsc(6+6*i,tn)*2*eta(i+1)/deltat+b3c(6+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))

epsnc(3+6*i,tn+1)=(epsnc(3+6*i,tn)*2*eta(i+1)/deltat+b3nc(3+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
epsnc(4+6*i,tn+1)=(epsnc(4+6*i,tn)*2*eta(i+1)/deltat+b3nc(4+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
epsnc(6+6*i,tn+1)=(epsnc(6+6*i,tn)*2*eta(i+1)/deltat+b3nc(6+6*i))/(2*(mu(i+1)+eta(i+1)/deltat))
if ((tn>=2).and.(tn<pocetkroku)) then
p(i+1,tn)=((eps(3+6*i,tn+1)-eps(3+6*i,tn-1))**2+(eps(4+6*i,tn+1)-eps(4+6*i,tn-1))**2+(eps(6+6*i,tn+1)-eps(6+6*i,tn-1))**2+&
2*((epsnc(3+6*i,tn+1)-epsnc(3+6*i,tn-1))**2+(epsnc(4+6*i,tn+1)-epsnc(4+6*i,tn-1))**2+(epsnc(6+6*i,tn+1)-epsnc(6+6*i,tn-1))**2)+&
2*((epsc(3+6*i,tn+1)-epsc(3+6*i,tn-1))**2+(epsc(4+6*i,tn+1)-epsc(4+6*i,tn-1))**2+(epsc(6+6*i,tn+1)-epsc(6+6*i,tn-1))**2))/&
deltat**2/4+&
p(i+1,tn)
endif
enddo
call vypocitejur(pocetvrstev,l,b3,ur)
call vypocitejur(pocetvrstev,l,b3nc,ur22)

!radial displacement calculation
theta = pi/2
urmemory = (ur(1)+ur(2))/2*1/4*sqrt(5.0d0/pi)*(3*cos(theta)*cos(theta)-1)+&
(ur22(1)+ur22(2))*1/4*sqrt(15.0d0/2/pi)*sin(theta)*sin(theta)

if (tn>pocetkroku-krokuvperiode) then
if (abs(urmemory) > urmaxmemory) then
urmaxmemory = abs(urmemory)
endif
if ((ur(1)+ur(2))/2<urmax) then
    urmax=(ur(1)+ur(2))/2
    tnmax=1.0d0*(tn-(pocetkroku-krokuvperiode))/krokuvperiode
endif
endif

if (tnmax.eq.1) then
tnmax = 0
endif
write(3,*), tn, (ur(1)+ur(2))/2  
end do
close(3)
write(5,*), eta(1), urmax
write(12,*), eta(1), tnmax
write(18,*), eta(1), urmaxmemory
!reversed eta in comparison to maxwell - need to adjust it for procedure
do i=1,pocetvrstev
eta(i)=1/eta(i)/4
enddo
call vypocitejpcelk(pcelk,pocetvrstev,pocetkroku,p,deltat,eta,r,krokuvperiode)
do i=1,pocetvrstev
eta(i)=1/eta(i)/4
enddo
print*, pcelk
write(2,*),eta(1),pcelk
!konec kelvin voigtova reologickeho modelu
!----------------------------------------------------------------------------------------------------------
enddo

close(1)
close(2)
close(4)
close(5)
close(11)
close(12)
close(17)
close(18)

deallocate(bncmemory)
deallocate(bcmemory)
deallocate(b2nc)
deallocate(b2c)
deallocate(a2)
deallocate(bmemory)
deallocate(b2)

deallocate(epsnc)
deallocate(epsc)
deallocate(epsmemory)
deallocate(b3nc)
deallocate(b3c)
deallocate(eps)
deallocate(a3)
deallocate(bmemory3)
deallocate(b3)
deallocate(p)

deallocate(ur22)
deallocate(ur)
deallocate(b)
deallocate(indx)
deallocate(al)
deallocate(a)

deallocate (eta)
deallocate (rtilded)
deallocate (rhatd)
deallocate (mu)
deallocate (rtilde)
deallocate (rhat)
deallocate (r)

end program naplnenimatice

!Numerical recipes subroutines
SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
INTEGER m1,m2,mp,mpl,n,np,indx(n)
REAL a(np,mp),al(np,mpl),b(n)
INTEGER i,k,l,mm
REAL dum
mm=m1+m2+1
!if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in banbks'
l=m1
do k=1,n
i=indx(k)
if(i.ne.k)then
dum=b(k)
b(k)=b(i)
b(i)=dum
endif
if(l.lt.n)l=l+1
do i=k+1,l
b(i)=b(i)-al(k,i-k)*b(k)
enddo
enddo
l=1
do i=n,1,-1
dum=b(i)
do k=2,l
dum=dum-a(i,k)*b(k+i-1)
enddo
b(i)=dum/a(i,1)
if(l.lt.mm) l=l+1
enddo
return
END subroutine
SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
INTEGER m1,m2,mp,mpl,n,np,indx(n)
REAL d,a(np,mp),al(np,mpl),TINY
PARAMETER (TINY=1.e-20)
INTEGER i,j,k,l,mm
REAL dum
mm=m1+m2+1
!if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
l=m1
do i=1,m1
do j=m1+2-i,mm
a(i,j-l)=a(i,j)
enddo
l=l-1
do j=mm-l,mm
a(i,j)=0.
enddo

enddo

d=1.
indx(k)=i
if(dum.eq.0.) a(k,1)=TINY
if(i.ne.k)then
d=-d
do j=1,mm
dum=a(k,j)
a(k,j)=a(i,j)
a(i,j)=dum
enddo
endif
do i=k+1,l
dum=a(i,1)/a(k,1)
al(k,i-k)=dum
do j=2,mm
a(i,j-1)=a(i,j)-dum*a(k,j)
enddo
a(i,mm)=0.
enddo
!enddo
return
END subroutine
subroutine bandecdp(a,n,m1,m2,np,mp,al,mpl,indx,d) !doubleprecision bandec
      INTEGER m1,m2,mp,mpl,n,np,indx(*)
      REAL*8 d,a(np,mp),al(np,mpl),TINY
      PARAMETER (TINY=0.d-20)
      INTEGER i,j,k,l,mm
      REAL*8 dum

      mm=m1+m2+1
      !if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
      l=m1
      do 13 i=1,m1
        do 11 j=m1+2-i,mm
          a(i,j-l)=a(i,j)
11      continue
        l=l-1
        do 12 j=mm-l,mm
          a(i,j)=0d0
12      continue
13    continue
      d=1d0
      l=m1
      do 18 k=1,n
        dum=a(k,1)
        i=k
        if(l.lt.n)l=l+1
        do 14 j=k+1,l
          if(abs(a(j,1)).gt.abs(dum))then
            dum=a(j,1)
            i=j
          endif
14      continue
        indx(k)=i
        if(dum.eq.0d0) a(k,1)=TINY
        if(i.ne.k)then
          d=-d
          do 15 j=1,mm
            dum=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=dum
15        continue
        endif
        do 17 i=k+1,l
          dum=a(i,1)/a(k,1)
          al(k,i-k)=dum
          do 16 j=2,mm
            a(i,j-1)=a(i,j)-dum*a(k,j)
16        continue
          a(i,mm)=0d0
17      continue
18    continue
      return
      END subroutine

      SUBROUTINE banbksdp(a,n,m1,m2,np,mp,al,mpl,indx,b)
      INTEGER m1,m2,mp,mpl,n,np,indx(*)
      REAL*8 a(np,mp),al(np,mpl),b(*)
      INTEGER i,k,l,mm
      REAL*8 dum

      mm=m1+m2+1
      !if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in banbks'
      l=m1
      do 12 k=1,n
        i=indx(k)
        if(i.ne.k)then
          dum=b(k)
          b(k)=b(i)
          b(i)=dum
        endif
        if(l.lt.n)l=l+1
        do 11 i=k+1,l
          b(i)=b(i)-al(k,i-k)*b(k)
11      continue
12    continue
      l=1
      do 14 i=n,1,-1
        dum=b(i)
        do 13 k=2,l
          dum=dum-a(i,k)*b(k+i-1)
13      continue
        b(i)=dum/a(i,1)
        if(l.lt.mm) l=l+1
14    continue
      return
      END subroutine
