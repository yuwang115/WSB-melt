
datadir = "/g/data/jk72/yw0666/data/"
outdir = "./vtuoutputs"
meshdb = "../mesh2D_WSB_BM3_refined"

-- ## Inversion regularisation
-- ## Lambda=1.0e10


-- ## Min SSA Critical Thickness
Hmin = 15

-- ## Min threshold value for Ice thickness (Real)
MINH=15

-- ## levels in the vertical 30
MLEV=15  

-- ## controlling steady state iterations
IMIN=10
IMAX=300

Tol=0.01
DPtol = 0.001

-- ## for block preconditioner 
blocktol=0.001

-- ##
name = "WSB_S0"

-- ## Mesh refinement parameters

-- ##  the name of the mesh to optimize
MESH_IN="mesh2D_WSB_BM3"

MESH_OUT="mesh2D_WSB_BM3_refined"

-- ## Tolerated errors on U and H
U_err= 6  -- 2
H_err= 500    -- 35.0


-- ## mesh size limits in different regions

-- ## absolute minimum mesh size
Mminfine=600.0

-- ## minimum mesh size far from grounding line (may be higher than
-- ## Mminfine to prevent detailed refinement in slow flowing far 
-- ## inland regions).  Set equal to Mminfine if you don't want this.
Mmincoarse=1000.0

-- ## Maximum mesh size far from the grounding line 15000
Mmaxfar=8000.0

-- ## Maximum mesh size close to the grounding line 2000
Mmaxclose=1000.0

-- ## maximum mesh size for ice shelves (set this to larger than
-- ## Mmaxfar if you want it to be ignored)
Mmaxshelf=1000.0

-- ## reference velocity used in refinement to set upper limit on
-- ## element size (together with distance from GL).  Sections of
-- ## grounding line with flow speeds higher than this limit will
-- ## have max mesh size Mmaxclose.  Max mesh is allowed to be
-- ## coarser for sections of slower flowing GL.  Set this very
-- ## small (e.g. 0.1) if you want it to be ignored.
refvel = 500 -- 700.0

-- ## The distance from grounding line at which the upper limit for
-- ## the maximum mesh size is reached
GLdistlim=300000.0

-- ## The distance from the boundary at which the upper limit for
-- ## the maximum mesh size is reached 80000
Bdistlim=70000.0

-- ## For distances beyond distlim, the minimum mesh size may be
-- ## increased towards Mmincoarse on this distance scale 400000
distscale=300000.0

-- TOPOGRAPHY_DATA="/g/data/jk72/cxz581/data/antarctic/bedmachineantarctica_v3_extended.nc"
TOPOGRAPHY_DATA="/g/data/jk72/cxz581/data/antarctic/bedmachineantarctica_v3_extended_watercolumn.nc"
VELOCITY_DATA="/g/data/jk72/cxz581/data/antarctic/antarctic_ice_vel_phase_map_v01_slim.nc"
DHDT_DATA = ""
SLIP_DATA = "/g/data/jk72/cxz581/data/antarctic/aa_v3_e8_l11_beta.nc"
SMB_DATA = "/g/data/jk72/cxz581/data/antarctic/smbref_1995_2014_mar.nc"
VISCOSITY_DATA = "/g/data/jk72/cxz581/data/antarctic/ant08_b2_future25-06_hist0001_mumean_pat94.nc"
-- ## BMB_DATA_SUSHEEL= "/g/data/jk72/cxz581/data/antarctic/bmb_susheel.nc"


-- ## Regularized the grounded mask, converting non-integers (decimals) present in the grounded mask to -1, 0 and 1.
function GroundedMaskRegularization(groundedmask)

      if (groundedmask > 0.5) then
        groundedmask = 1
      elseif (groundedmask < -0.5) then
        groundedmask = -1
      else
      groundedmask = 0
      end
    return groundedmask
end 


function IfThenElse(condition,t,f) 
  if condition then
    return t
  else
    return f
  end 
end


function BedDeepening(groundedmask, bedrock)
  if groundedmask == -1 then
    bedrock = bedrock - 10
  end
  return bedrock
end



RandomNumbers = {298,297,304,291,292,288,305,299,296,299,297,299,297,301,296,306,290,288,288,287,294,295,293,302,299,302,305,306,290,288,300,287,297,297,304,296,294,300,301,296,293,289,298,291,286,301,291,295,300,293,301,294,300,300,295,286,292,294,291,290,303,295,304,294,302,294,302,301,293,290,302,305,292,300,295,303,302,289,304,306,296,304,298,289,290,294,301,303,302,292,297,287,288,288,300,296,289,296,289,287,303,297,305,300,298,303,304,306,286,304,298,306,297,296,302,290,296,304,298,303,301,298,291,299,287,299,299,301,304,306,302,298,305,298,286,288,304,296,303,290,297,299,286,298,293,287,296,290,288,290,289,289,286,299,291,297,300,296,297,295,288,296,303,304,291,290,297,299,294,290,305,287,288,288,289,299,298,287,305,301,301,287,304,305,306,304,302,296,289,294,288,286,305,292,292,292,295,299,286,303,297,303,293,295,287,289,299,292,304,288,306,297,300,306,292,294,295,302,303,288,289,293,287,296,293,289,290,305,300,295,305,288,301,301,297,289,298,292,288,290,304,287,291,287,295,286,304,290,287,292,295,288,306,292,292,287,292,286,296,301,299,287,287,302,305,297,288,303,293,292,301,286,287,300,298,297,301,300,302,292,300,297,294,287,302,293,298,301,288,288,297,296,304,302,301,287,287,287,302,305,300,288,301,288,288,299,292,299,301,298,301,290,301,306,304,287,293,293,300,298,302,293,290,287,302,290,294,297,290,299,296,289,302,288,292,290,297,287,294,288,288,302,292,298,306,295,300,301,295,299,288,305,289,291,302,296,302,294,291,286,300,295,295,298,287,292,302,300,288,288,287,286,294,299,301,297,288,299,288,288,288,288,289,290,292,292,290,291,304,300,297,289,290,287,305,300,297,292,289,299,306,289,291,294,287,300,294,306,294,299,289,294,289,301,304,293,300,292,297,303,298,293,292,295,294,293,297,301,294,295,288,286,292,292,299,306,305,295,291,302,301,301,301,288,300,295,290,288,303,289,289,299,304,296,300,289,306,297,300,286,302,301,288,297,292,297,294,294,289,291,286,305,299,305,289,305,302,298,295,291,301,290,287,302,300,301,299,294,294,303,292,303,302,303,296,299,305,295,287,304}
RandomNumberIndex = 1
function calculateTimePoint(Time, StartFrom)
  local TimePoint = Time + StartFrom
  if TimePoint <= 306 then
      return TimePoint
  else
      local randomNumber = RandomNumbers[RandomNumberIndex]
      RandomNumberIndex = RandomNumberIndex + 1
      return randomNumber
  end
end