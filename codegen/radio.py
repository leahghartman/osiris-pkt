# Set up the combination of field components, with the E fields listed first
import itertools
#do exclusive E_field combs
combs1=[]
for i in range(1,4):
    combs1=combs1+list(itertools.combinations([11,12,13], i))
#do E_field+B_field combs
combs2=[]
for i in range(1,7):
    combs2=combs2+list(itertools.combinations([11,12,13,21,22,23], i))
#remove exclusive E_field combs from E_field+B_field combs
for a in combs1:
    if a in combs2: combs2.remove(a)
#Now exclusive efield combs come first
combs=combs1+combs2

# Open destination file for macros
file= open("os-macros.f03","w")
body = ""
nc = "#define __NC__ "
first_b = True

# Write out the macros for each combination
for combi in combs:

    #auxilliary variables
    funcname="#define FNAME(a) a ## _"
    ind = ""
    count=1+0
    flag=True

    #loop through the components of the current combination
    for i in combi:
        #get the spatial component
        j=int(str(i)[1])
        #distinguish between the Electic and the Magnetic field
        if(int(str(i)[0])==1):
            fld="E"
        else:
            if first_b:
                body = body + "#ifdef __HAS_RAD_BFLD__\n\n"
                first_b = False
            fld="B"
        ind = ind + "#define __I" + fld + str(i)[1] + "__ " + str(count) + "\n"
        #keep track of the current component
        count=count+1
        #update the subroutine name
        funcname=funcname+fld+str(i)[1]

    body = body + funcname + "\n"
    body = body + nc + str(count - 1) + "\n"
    body = body + ind
    body = body + "#include __FILE__\n\n"

body = body + "#endif"

file.write(body)
file.close()


# Now generate the selection functions for filling emf
invars="  real(p_k_part), dimension(3), intent(in) :: x, xprev\n"
invars=invars+"  real(p_k_part), dimension(p_p_dim), intent(in) :: beta, betaprev\n"
invars=invars+"  real(p_k_part), intent(in) :: charge\n"
invars=invars+"  real(p_double), intent(in) :: time\n"
invars=invars+"  integer, intent(in) :: tid\n\n  ! Do nothing\n\n"
dash = "!-----------------------------------------------------------------------------------------\n"

def manyifs(decomp):
    inbod=dash+"subroutine "+"comp_select_"+decomp+"_det( this )\n\n"
    inbod=inbod+"  implicit none\n\n"
    inbod=inbod+"  class("+decomp+"_Detector), intent(inout) :: this\n\n"
    inbod=inbod+"  integer :: sele\n\n"
    inbod=inbod+"  sele = 0\n"
    inbod=inbod+"  if (this%cmpts(6)) sele = sele + 32\n"
    inbod=inbod+"  if (this%cmpts(5)) sele = sele + 16\n"
    inbod=inbod+"  if (this%cmpts(4)) sele = sele + 8\n"
    inbod=inbod+"  if (this%cmpts(3)) sele = sele + 4\n"
    inbod=inbod+"  if (this%cmpts(2)) sele = sele + 2\n"
    inbod=inbod+"  if (this%cmpts(1)) sele = sele + 1\n\n"

    if decomp=="cart":
        sp=""
        nn=1
        apprx=["approx_"]
    else:
        inbod=inbod+"  if (this%approx) then\n\n"
        sp="  "
        nn=2
        apprx=["approx_",""]

    for n in range(nn):
        inbod=inbod+sp+"  select case(sele)\n"
        for i in range(1,64):
            if (i==8): inbod=inbod+"#ifdef __HAS_RAD_BFLD__\n"
            ifstat=sp+"  case({})\n".format(i)+sp+"    this % fill_emf_loc => "
            name="calc_comp_"+apprx[n]
            if(1&i):
                name=name+"E1"
            if(2&i):
                name=name+"E2"
            if(4&i):
                name=name+"E3"
            if(8&i):
                name=name+"B1"
            if(16&i):
                name=name+"B2"
            if(32&i):
                name=name+"B3"
            inbod=inbod+ifstat+name+"\n"
        inbod=inbod+"#endif\n"

        inbod=inbod+sp+"  case default\n"+sp+"    this % fill_emf_loc => calc_dflt\n"
        inbod=inbod+sp+"  end select\n\n"
        if decomp=="sphe":
            if n==0:
                inbod=inbod+"  else\n\n"
            else:
                inbod=inbod+"  endif\n\n"

    inbod=inbod+"end subroutine "+"comp_select_"+decomp+"_det\n"+dash
    return inbod


pre=dash
pre=pre+"subroutine calc_dflt(this,x,beta,time,xprev,betaprev,charge,tid)\n"
pre=pre+"  implicit none\n"
pre=pre+"  class(sphe_Detector), intent(inout) :: this\n"
pre=pre+invars
pre=pre+"end subroutine calc_dflt\n"
pre=pre+dash+"\n"

mainbod=pre
mainbod=mainbod+manyifs("sphe")+"\n"

file= open("os-sphe-sel.f03","w")
file.write(mainbod)
file.close()



mainbod=pre
mainbod=mainbod+manyifs("cart")+"\n"

file= open("os-cart-sel.f03","w")
file.write(mainbod.replace("sphe","cart"))
file.close()
