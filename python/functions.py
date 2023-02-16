import maya.cmds
import maya.mel as mel
import das
from itertools import cycle


"""
curve = maya.cmds.ls(sl=True)[0]
params = GetRange(count=16)
joints = MakeJoints(name = "rail_crv_OB", count = len(params), intm = False, addAtts=None)
obs = maya.cmds.ls(sl=True)
aim = maya.cmds.ls(sl=True)[0]
surf = maya.cmds.ls(sl=True)[0]

pociNodes = CurveParamsToTranslate(params=params, curve=curve, top=True, obs=obs)
rNodes = GetRotationFromPoci(pociNodes=pociNodes, nodes=obs, aim=aim)
rots = rNodes["comp"]

ConnectAtts(rots, "outputRotate", obs, "r", verbose=False)

uvVal = []
yuu = 2
vee = 17
for u in range(1,yuu):
        for v in range(1,vee):
            uv = [u / float(yuu), v / float(vee)]
            uvVal.append(uv)    

cfsiNodes = GetCurveFromSurfaceNodes(uvVal, surf, 0)
infoNodes = GetCurveInfo(cfsiNodes)  
baseLen = []
for i in range(len(infoNodes)):   
    baseLen.append(1.0 / maya.cmds.getAttr("%s.arcLength" % infoNodes[i]))
multNodes = GetMulti(cfsiNodes, baseLen)
ConnectAtts(infoNodes, "arcLength", multNodes, "input1")
ConnectAtts(multNodes, "output", obs, "sz")
ConnectAtts(multNodes, "output", obs, "sy")
blendos = maya.cmds.ls("blend_*_mat")


for x in blendos:
    connex = maya.cmds.listConnections("%s.wtMatrix[1].matrixIn" % x, s=True, p=True)[0]
    multi = maya.cmds.createNode("multMatrix", n=x.replace("blend", "multi"))
    #pMulti = maya.cmds.createNode("multMatrix", n=x.replace("blend", "inverseMulti"))
    maya.cmds.connectAttr(connex, "%s.matrixIn[0]" % multi)
    #maya.cmds.connectAttr("locator1.worldMatrix[0]", "%s.matrixIn[0]" % pMulti)
    maya.cmds.connectAttr("locator1.worldMatrix[0]", "%s.matrixIn[1]" % multi)
    #maya.cmds.connectAttr("locator1.parentInverseMatrix[0]", "%s.matrixIn[1]" % pMulti)
    #maya.cmds.connectAttr("%s.matrixSum" % pMulti, "%s.wtMatrix[1].matrixIn" % x, f=True)
    maya.cmds.connectAttr("%s.matrixSum" % multi, "%s.wtMatrix[1].matrixIn" % x, f=True)

for x in blendos:
    connex = maya.cmds.listConnections("%s.wtMatrix[1].matrixIn" % x, s=True, p=True)[0]
    pMulti = maya.cmds.createNode("multMatrix", n=x.replace("blend", "inverseMulti"))
    maya.cmds.connectAttr(connex, "%s.matrixIn[0]" % multi)
    maya.cmds.connectAttr("wrap_ctl_gp.parentInverseMatrix[0]", "%s.matrixIn[1]" % pMulti)
    maya.cmds.connectAttr("wrap_ctl_gp.inverseMatrix[0]", "%s.matrixIn[2]" % pMulti)
    maya.cmds.connectAttr("wrap_ctl.worldMatrix[0]", "%s.matrixIn[3]" % pMulti)
    maya.cmds.connectAttr("%s.matrixSum" % pMulti, "%s.wtMatrix[1].matrixIn" % x, f=True)
    
"""


import math


def MakeCtl(ob, name=None, meshDisplay=None, vis=None, guts=None):
    # A BF style control object
    name = ob if name is None else name
    size = maya.cmds.getAttr("%s.boundingBoxSize" % ob)[0]
    rad = math.sqrt(size[0]**2 + size[2]**2) / 2
    crcle, makeCircle = maya.cmds.circle(name="%s_ctl" % ob, r=rad, nry=1, nrx=0, nrz=0)
    gp = maya.cmds.group(crcle, name="%s_ctl_gp" % ob)
    pcon = maya.cmds.parentConstraint(ob, gp, mo=0)
    maya.cmds.delete(pcon)
    mat = maya.cmds.createNode("decomposeMatrix")
    mult = maya.cmds.createNode("multMatrix")
    maya.cmds.connectAttr("%s.worldMatrix[0]" % crcle, "%s.matrixIn[0]" % mult)
    maya.cmds.connectAttr("%s.parentInverseMatrix[0]" % crcle, "%s.matrixIn[1]" % mult)
    maya.cmds.connectAttr("%s.matrixSum" % mult, "%s.inputMatrix" % mat)
    maya.cmds.connectAttr("%s.outputTranslate" % mat, "%s.t" % ob)
    maya.cmds.connectAttr("%s.outputRotate" % mat, "%s.r" % ob)
    maya.cmds.connectAttr("%s.outputScale" % mat, "%s.s" % ob)
    if meshDisplay:
        maya.cmds.addAttr(crcle, ln="meshDisplay", k=True,  at="long", min=0, max=2, dv=2 )
        maya.cmds.connectAttr("%s.meshDisplay" % crcle, "%s.overrideEnabled" % ob)
        maya.cmds.connectAttr("%s.meshDisplay" % crcle, "%s.overrideDisplayType" % ob)
    if vis:
        if maya.cmds.objExists(vis):
            maya.cmds.addAttr(crcle, ln="vis", k=True,  at="long", min=0, max=1, dv=1 )
            maya.cmds.connectAttr("%s.vis" % crcle, "%s.v" % vis)
    if guts:
        if maya.cmds.objExists(guts):
            maya.cmds.addAttr(crcle, ln="guts", k=True,  at="long", min=0, max=1, dv=1 )
            maya.cmds.connectAttr("%s.guts" % crcle, "%s.v" % guts)

   
def SuperGroupIt(obs):
    # Group based on object name
    gp = ""
    if len(obs) == 0:
        gp = mel.eval("CreateEmptyGroup")
    else:
        if isinstance(obs, str):
            obs = [obs]
        groupName = "%s_gp" % obs[0]
        gp = maya.cmds.group(obs, name=groupName)
    return gp


def GetRange(lo=0.0, hi=1.0, count=100, precision = 3):
    params = []
    normal = float(hi - lo)
    for i in range(count):
        result = round(normal / (count-1) * i + float(lo), precision)
        params.append(result)
    return params


def Get2dRange(yuu=2, vee=8):
    range2d = []
    for u in range(1,yuu):
        for v in range(1,vee):
            uv = [u / float(yuu), v / float(vee)]
            range2d.append(uv)     
    return range2d


def MakeNewCurve(name="null_func_lr", fromSelection=True, points=[], degree = 3):
    # takes a single open edge selection and returns a curve and vertex worldspace positions
    pos = []
    split = name.rsplit("_", 1)
    name = "%s_crv_%s" % (split[0], split[1])

    if fromSelection:
        selVer = maya.cmds.polyListComponentConversion(fe=True, tv=True)
        verts = maya.cmds.ls(selVer, fl=True)
        for v in verts:
            pos.append(maya.cmds.pointPosition(v))
        newCurve = maya.cmds.polyToCurve(form=0, degree=degree, ch=False, n=name)[0]

    else:
        pos = points
        if not pos:
            return None
        newCurve = maya.cmds.curve(p=pos, d=degree, n=name)
            
    return newCurve, pos


def MakeSampleCurve(inputCurve = "null", spans = 4):
    #takes an input curve and returns a simplified curve with JNT name
    split = inputCurve.rsplit("_", 2)
    name = "%s_jnt_crv_%s" % (split[0], split[2])
    sampleCurve = maya.cmds.rebuildCurve(inputCurve, rpo = 0, ch = 0, d = 3, s = spans, n = name)[0]
    return sampleCurve


def MakeJoints(name = "null_obj_OB", count = 1, intm = True, addAtts=None):
    data = das.Struct()
    data.groups = []
    data.joints = []
    for num in range(count):
        val = ""
        if count > 1:
            val = "_%02d" % num 
        jntName = "%s%s_jnt" % (name, val)
        topName = "%s%s_grp" % (name, val)
        topGp = maya.cmds.group(n=topName, em=True)
        parent = topGp
        if intm:
            midName = "%s%s_int" % (name, val)
            parent = maya.cmds.group(p=topGp, n=midName, em=True) 
        newJoint = maya.cmds.createNode("joint", p=parent, n=jntName) 
        if addAtts:
            for attr in addAtts:
                maya.cmds.addAttr(newJoint, ln=attr, at="double", k=True)
        data.joints.append(newJoint) 
        data.groups.append(topGp) 
    return data


def GetCurvePoints(curve=None, ch=False):
    split = curve.rsplit("_", 2) or curve
    data = das.Struct()

    data.points = []
    data.curve = maya.cmds.listRelatives(curve, c = True, ni = True, type = "nurbsCurve")[0]
    data.spans = maya.cmds.getAttr("%s.spans" % data.curve) + maya.cmds.getAttr("%s.degree" % data.curve)
    data.node = maya.cmds.createNode("curveInfo", ss= True, n="%s_info" % split[0])

    maya.cmds.connectAttr("%s.worldSpace[0]" % data.curve, "%s.inputCurve" % data.node)
    for i in range(data.spans):
        pt = maya.cmds.getAttr("%s.controlPoints[%s]" % (data.node, i))[0]
        data.points.append(pt)

    if not ch: 
        maya.cmds.delete(data.node)
        data.node = None

    return data


def ConnectCurvePointsToTranslate(curveInfo=None, nodes=[]):
    for i in range(len(nodes)):
        maya.cmds.connectAttr("%s.controlPoints[%s]" % (curveInfo, i), "%s.t" % nodes[i], f=True)


def CurveParamsToTranslate(params=[0.0,0.5,1.0], curve=None, top=False, obs=None):
    #Attach objects to curve params or return poci Nodes
    pociNodes = []
    split = curve.rsplit("_", 2) or curve
    curv = maya.cmds.listRelatives(curve, c = True, ni = True, type = "nurbsCurve")[0]
    for i in range(len(params)):
        poci = maya.cmds.createNode("pointOnCurveInfo", ss = True, n = "%s_%02d_poci" % (split[0], i))
        maya.cmds.setAttr("%s.parameter" % poci, params[i])
        maya.cmds.setAttr("%s.turnOnPercentage" % poci, top)
        maya.cmds.connectAttr("%s.worldSpace[0]" % curv, "%s.inputCurve" % poci)
        pociNodes.append(poci)
    if obs:
        if len(obs) != len(pociNodes):
            print ("Objects and pointOnCurve nodes are not the same count")
        try:
            AttachObsToInfoNodes(obs, pociNodes, t=False, r=False)
        except:
            print("Well that didn't work")
    return pociNodes


def GetRotationFromPoci(pociNodes=None, nodes=None, aim=None):
    data = das.Struct()
    data.fbfm = []
    data.dcmp = []
    data.vpro = []
    data.mult = []
    data.comp = []

    for i in range(len(pociNodes)):
        fbfm = maya.cmds.createNode("fourByFourMatrix", ss=True, n=pociNodes[i].replace("poci", "fbfm"))
        dcmp = maya.cmds.createNode("decomposeMatrix",  ss=True, n=pociNodes[i].replace("poci", "dcmp"))
        mult = maya.cmds.createNode("multMatrix",       ss=True, n=pociNodes[i].replace("poci", "mult"))
        vpro = maya.cmds.createNode("vectorProduct",    ss=True, n=pociNodes[i].replace("poci", "vpro"))
        comp = maya.cmds.createNode("decomposeMatrix",  ss=True, n=pociNodes[i].replace("poci", "comp"))
        maya.cmds.setAttr("%s.operation" % vpro, 2)
        data.vpro.append(vpro)
        data.fbfm.append(fbfm)
        data.dcmp.append(dcmp)
        data.mult.append(mult)
        data.comp.append(comp)

    ConnectAtts(pociNodes,      "tangent",                  data.vpro,      "input1")
    ConnectAtts(data.dcmp,      "outputTranslate",          data.vpro,      "input2")
    ConnectAtts(nodes,          "worldInverseMatrix[0]",    data.mult,      "matrixIn[0]")
    ConnectAtts(aim,            "worldMatrix[0]",           data.mult,      "matrixIn[1]")
    ConnectAtts(data.mult,      "matrixSum",                data.dcmp,      "inputMatrix")
    ConnectAtts(pociNodes,      "tangentX",                 data.fbfm,      "in00")
    ConnectAtts(pociNodes,      "tangentY",                 data.fbfm,      "in01")
    ConnectAtts(pociNodes,      "tangentZ",                 data.fbfm,      "in02")
    ConnectAtts(data.vpro,      "outputX",                  data.fbfm,      "in10")
    ConnectAtts(data.vpro,      "outputY",                  data.fbfm,      "in11")
    ConnectAtts(data.vpro,      "outputZ",                  data.fbfm,      "in12")
    ConnectAtts(data.dcmp,      "outputTranslateX",         data.fbfm,      "in20")
    ConnectAtts(data.dcmp,      "outputTranslateY",         data.fbfm,      "in21")
    ConnectAtts(data.dcmp,      "outputTranslateZ",         data.fbfm,      "in22")
    ConnectAtts(data.fbfm,      "output",                   data.comp,      "inputMatrix")

    return data


def GetClosestPointsOnCurve(points, sampleCurve):
    params = []
    curv = maya.cmds.listRelatives(sampleCurve, c = True, ni = True, type = "nurbsCurve")[0]
    npoc = maya.cmds.createNode("nearestPointOnCurve", ss = True)
    maya.cmds.connectAttr("%s.worldSpace[0]" % curv, "%s.inputCurve" % npoc)
    for one in points:
        maya.cmds.setAttr("%s.inPosition" % npoc, one[0], one[1], one[2])
        value = maya.cmds.getAttr("%s.result.parameter" % npoc )
        params.append(value)
    maya.cmds.delete(npoc)
    return params


def GetCurveInfoNodes(cfsiNodes):
    # does it have a worldspace or an outputcurve?
    # This only accepts curveFromSurface nodes, should also accept curves
    infoNodes =[]
    output = "outputCurve"
    suffix = cfsiNodes[0].split("_")[-1]
    for i in range(len(cfsiNodes)):
        info = maya.cmds.createNode("curveInfo", ss = True, name = cfsiNodes[i].replace(suffix, "info")) 
        maya.cmds.connectAttr("%s.%s" % (cfsiNodes[i], output), "%s.inputCurve" % info, f = True)
        infoNodes.append(info)
    return infoNodes


def GetMulti(nodes, defaultVal):
    multNodes = []
    suffix = nodes[0].split("_")[-1]
    for i in range(len(nodes)):
        mult = maya.cmds.createNode("multDoubleLinear", ss = True, name = nodes[i].replace(suffix, ""))
        maya.cmds.setAttr("%s.input2" % mult, defaultVal[i])
        multNodes.append(mult)
    return multNodes


def GetCurveFromSurfaceNodes(params, surface, uv, relative=False):
    cfsINodes = []
    surf = maya.cmds.listRelatives(surface, c=True, ni=True, type="nurbsSurface")[0] 
    for i in range(len(params)):
        value = (uv *-1 )+ 1
        parameter = params[i][value]
        cfsI = maya.cmds.createNode("curveFromSurfaceIso", ss = True, n = "%s_%02d_cfsInfo" % (surface, i))
        maya.cmds.connectAttr("%s.worldSpace[0]" % surf, "%s.inputSurface" % cfsI)
        maya.cmds.setAttr("%s.isoparmDirection" % cfsI, uv)
        maya.cmds.setAttr("%s.relativeValue" % cfsI, relative)
        maya.cmds.setAttr("%s.isoparmValue" % cfsI, parameter)
        cfsINodes.append(cfsI)
    return cfsINodes
   

def GetClosestPointsOnSurface(points, sampleSurface):
    cposNodes = []
    split = sampleSurface.split("_")
    surf = maya.cmds.listRelatives(sampleSurface, c=True, ni=True, type="nurbsSurface")[0]  
    for i in range(len(points)):  
        cPos = maya.cmds.createNode("closestPointOnSurface", ss=True, n="%s_%02d_cPos" % (split[0], i))
        maya.cmds.connectAttr("%s.worldSpace[0]" % surf, "%s.inputSurface" % cPos)
        maya.cmds.setAttr("%s.inPosition" % cPos, points[i][0], points[i][1], points[i][2])
        u = maya.cmds.getAttr("%s.result.parameterU" % cPos )
        v = maya.cmds.getAttr("%s.result.parameterV" % cPos )
        cposNodes.append(cPos)
    return cposNodes


# def GetClosestPointsOnMesh(points, sampleSurface):
#     cposNodes = []
#     split = sampleSurface.split("_")
#     surf = maya.cmds.listRelatives(sampleSurface, c=True, ni=True, type="nurbsSurface")[0]  
#     for i in range(len(points)):  
#         cPos = maya.cmds.createNode("closestPointOnSurface", ss=True, n="%s_%02d_cPos" % (split[0], i))
#         maya.cmds.connectAttr("%s.worldSpace[0]" % surf, "%s.inputSurface" % cPos)
#         maya.cmds.setAttr("%s.inPosition" % cPos, points[i][0], points[i][1], points[i][2])
#         u = maya.cmds.getAttr("%s.result.parameterU" % cPos )
#         v = maya.cmds.getAttr("%s.result.parameterV" % cPos )
#         cposNodes.append(cPos)
#     return cposNodes


def GetPointOnSurfaceInfo(params, sampleSurface, percentage=True):
    posiNodes = []
    surf = maya.cmds.listRelatives(sampleSurface, c=True, ni=True, type="nurbsSurface")[0] 
    split = sampleSurface.rsplit("_", 1)
    for i in range(len(params)):
        posi = maya.cmds.createNode("pointOnSurfaceInfo", n = "%s_%02d_posi" % (split[0], i), ss= True)    
        maya.cmds.setAttr("%s.parameterU" % posi, params[i][0])
        maya.cmds.setAttr("%s.parameterV" % posi, params[i][1])
        maya.cmds.setAttr("%s.turnOnPercentage" % posi, percentage)
        maya.cmds.setAttr("%s.parameterU" % posi, k = True)
        maya.cmds.setAttr("%s.parameterV" % posi, k = True)
        maya.cmds.connectAttr("%s.worldSpace[0]" % surf, "%s.inputSurface" % posi) 
        posiNodes.append(posi)
    return posiNodes



def NormalAimConstraint(obs, points, up, aimVec=[1,0,0], upVec=[0,0,1], keepHistory=True):
    # Aim many objects to many points or many objects to one point
    nodes = []
    for i in range(len(obs)):
        if len(obs) == len(points):
            target = points[i]
        else:
            target = points
        aimCon = maya.cmds.aimConstraint(target, obs[i], mo=False, aim=aimVec, u=upVec, wut="object", wuo=up)
        if keepHistory:
            nodes.append(aimCon[0])
        else:
            maya.cmds.delete(aimCon)
    return nodes


def AlignObsToSurfaceParams(obs, posiNodes):
    #requires posi nodes and objects
    dcmpNodes = []
    fbfmNodes = []
    split = posiNodes[0].rsplit("_", 2)
    for i in range(len(posiNodes)): 
        if maya.cmds.objExists(obs[i]):
            posi = posiNodes[i]   
            dcmp = maya.cmds.createNode("decomposeMatrix", ss = True, name = "%s_dcmp%02d_%s" % (split[0], i, split[2])) 
            fbfm = maya.cmds.createNode("fourByFourMatrix", ss = True, name = "%s_fbfm%02d_%s" % (split[0], i, split[2])) 
            for input, output in  zip("012", "XYZ"):
                maya.cmds.connectAttr("%s.normal%s" % (posi, output), "%s.in0%s" % (fbfm, input))
                maya.cmds.connectAttr("%s.tangentU%s" % (posi, output.lower()), "%s.in1%s" % (fbfm, input))
                maya.cmds.connectAttr("%s.tangentV%s" % (posi, output.lower()), "%s.in2%s" % (fbfm, input))
            maya.cmds.connectAttr("%s.output" % fbfm, "%s.inputMatrix" % dcmp)
            maya.cmds.connectAttr("%s.outputRotate" % dcmp, "%s.r" % obs[i])
            dcmpNodes.append(dcmp)
            fbfmNodes.append(fbfm)
    return dcmpNodes, fbfmNodes


def ConnectSurfaceToScale(sampleSurface, uvVal, obs, axis, uvScale=0, relative=False):
    #Requires a sample surface, a list of obs, uvCoordinates and which sxes you wish to scale (accepts "x", "yz", "xyz" etc...)
    for xyz in axis:
        if xyz not in "xyz":
            print("scale only accepts x, y or z, got %s" % axis)
            return
    baseLen = []
    cfsiNodes = GetCurveFromSurfaceNodes(uvVal, sampleSurface, uvScale, relative)
    infoNodes = GetCurveInfoNodes(cfsiNodes)  
    for i in range(len(infoNodes)):   
        baseLen.append(1.0 / maya.cmds.getAttr("%s.arcLength" % infoNodes[i]))
    multNodes = GetMulti(cfsiNodes, baseLen)
    ConnectAtts(infoNodes, "arcLength", multNodes, "input1")
    for xyz in axis:
        ConnectAtts(multNodes, "output", obs, "s%s" % xyz)


def AttachObsToInfoNodes(obs, infoNodes, t=False, r=False):
    if t:
        ConnectAtts(infoNodes, "result.position", obs, "t")
    if r:
        dcmpNodes, fbfmNodes = AlignObsToSurfaceParams(obs, infoNodes)
        return (dcmpNodes, fbfmNodes)


def ConnectAtts(sourceNodes, sourceAtt, destNodes, destAtt=None, verbose=False):
    #connects one to many or n to n attributes, and makes attribute if it doesn"t exist
    single = False
    inn = sourceNodes
    out = destNodes
    destAtt = sourceAtt if destAtt is None else destAtt
    
    if not isinstance(destNodes, (list, tuple,)):
        out = [destNodes]

    if isinstance(sourceNodes, (list, tuple, set)):   
        sLen = len(sourceNodes)
        dLen = len(destNodes)
        if sLen != dLen:
            print("in Nodes and out nodes have different lengths.  (%s -> %s) Function works with one to many, or same number to same number" % (sLen, dLen))
            return
    else:
        single = True
        inn = [sourceNodes]

    node = inn[0]

    for i in range(len(out)):
        outNode = destNodes[i]
        if not single:
            node = inn[i] 
        try:
            maya.cmds.connectAttr("%s.%s" % (node, sourceAtt), "%s.%s" % (outNode, destAtt))
        except:
            print ("%s.%s could not connect to %s.%s" % (node, sourceAtt, outNode, destAtt))


def ConnectObsToCurve(obs, curve):
    # Connect the points of an existing curve to the translate of objects
    curv = maya.cmds.listRelatives(curve, c = True, ni = True, type = "nurbsCurve")[0] 
    for i in range(len(obs)):
        dcmp = maya.cmds.createNode("decomposeMatrix", ss = True, name = "%s_%02d_decomp" % (obs[i], i)) 
        maya.cmds.connectAttr("%s.worldMatrix[0]" % obs[i], "%s.inputMatrix" % dcmp)
        maya.cmds.connectAttr("%s.outputTranslate" % dcmp, "%s.controlPoints[%s]" % (curv, i)) 


def MatrixConstraint(obs):
    source = obs[0]
    for obj in obs[1:]:
        try:
            mM = maya.cmds.createNode("multMatrix", ss = True)
            decomp = maya.cmds.createNode("decomposeMatrix", ss = True)
            maya.cmds.connectAttr("%s.worldMatrix[0]" % source, "%s.matrixIn[0]" % mM)
            maya.cmds.connectAttr("%s.parentInverseMatrix[0]" % obj, "%s.matrixIn[1]" % mM)
            maya.cmds.connectAttr("%s.matrixSum" % mM, "%s.inputMatrix" % decomp)
            maya.cmds.connectAttr("%s.outputTranslate" % decomp, "%s.t" % obj)
            maya.cmds.connectAttr("%s.outputRotate" % decomp, "%s.r" % obj)
            maya.cmds.connectAttr("%s.outputScale" % decomp, "%s.s" % obj)
        except:
            print ("%s failed" % obj)


def CreateBlendShapes(names=[], targets=[]):   
    # Give it a target mesh and list of blendShape names to initialize blendshape setup
    targetShapes=[]
    blendShapes=[]

    for i in range(len(names)):
        new=maya.cmds.duplicate(targets[0], name="%s_%s" % (names[i], targets[0]))[0]
        targetShapes.append(new)
    for target in targets:   
        split=target.rsplit("_",2)  
        bs=maya.cmds.blendShape(target, foc=True, n="%s_BS" % (targets[0]))[0]
        for i in range(len(targetShapes)):
            maya.cmds.blendShape (bs, edit=True, foc=True, t=(target, i, targetShapes[i], 1)) 
        blendShapes.append(bs)
    return blendShapes


def WeldTheBitches():
    #Copy skin weights from first vert to other verts
    verts = maya.cmds.ls(sl=True, fl = True)
    maya.cmds.select(verts[0])
    mel.eval("artAttrSkinWeightCopy")
    maya.cmds.select(verts)
    mel.eval("artAttrSkinWeightPaste")
    

def CreateCurveOrientedLocators(sampleCurve, aimAt, upLoc):
    #add a handler to accept either a curve or an array of points
    split = sampleCurve.rsplit("_",2)
    points = GetCurvePoints(sampleCurve)
    newLocators = []
    newInts = []
    for i in range(len(points)): 
        pt = points[i]
        newLo = maya.cmds.spaceLocator(n = "%s%02d_%s_loc" % (split[0], i, split[2]))[0]
        gpName = "%sGP" % newLo
        intName = "%sINTGP" % newLo
        intGp = maya.cmds.group(n = intName)
        newGp = maya.cmds.group(n = gpName)
        maya.cmds.select(aimAt, newGp)
        SnipSnap()
        maya.cmds.move(pt[0],pt[1],pt[2], intGp, ws = True)
        newLocators.append(newLo)   
        newInts.append(intGp)    
    aimNodes = NormalAimConstraint(obs = newInts, points = aimAt, up = upLoc, dir = (-1,0,0), keepHistory = False)
    maya.cmds.delete(aimNodes)
    return newLocators


def SnipSnap(fromNodes = None, toNodes = None):
    # Temp parent constraints to move obs from one place to another
    if fromNodes == None:
        fromNodes = [maya.cmds.ls(sl=True)[0]]
    if toNodes == None:
        toNodes = [maya.cmds.ls(sl=True)[1]]
    for inc, out in zip(cycle(fromNodes), toNodes):
        trash = maya.cmds.parentConstraint(inc, out, mo = 0)
        maya.cmds.delete(trash)


def SuperConnectO():
    # connect attributes based on channel box selection
    obs = maya.cmds.ls(sl=True)
    selected = mel.eval("selectedChannelBoxAttributes")
    for one in selected:
        for i in range(1, len(obs)):
            maya.cmds.connectAttr("%s.%s" % (obs[0], one), "%s.%s" % (obs[i], one))


def CopyCurveShape(keepHistory=True):
    # Attach the shape of one curve to many curves.  Useful for controllers which only need one shape input
    obs = maya.cmds.ls(sl = True)
    source = maya.cmds.listRelatives(obs[0], c = True, ni = True, type = "nurbsCurve")[0] 
    for one in obs[1:]:
        curv = maya.cmds.listRelatives(one, c = True, ni = True, type = "nurbsCurve")[0]
        maya.cmds.connectAttr("%s.worldSpace[0]" % source, "%s.create" % curv, f = True)
        if not keepHistory:
            maya.cmds.disconnectAttr("%s.worldSpace[0]" % source, "%s.create" % curv)


def JointsOnSurface(yuu=2, vee=8, scale=None, uvScale=0):
    obs = maya.cmds.ls(sl=True)
    for sampleSurface in obs:

        fbfmNodes = []
        decompNodes = []
        uvVal = Get2dRange(yuu, vee)

        posiNodes = GetPointOnSurfaceInfo(uvVal, sampleSurface)
        data = MakeJoints(name=sampleSurface, count=len(uvVal), intm=False)
        groups = data.groups
        # joints = data.joints

        for i in range(len(posiNodes)):
            fbfmt = maya.cmds.createNode("fourByFourMatrix", ss = True, n = posiNodes[i].replace("posi", "fbfm"))
            decomp = maya.cmds.createNode("decomposeMatrix", ss= True, n = posiNodes[i].replace("fbfm", "dcmp"))
            fbfmNodes.append(fbfmt)
            decompNodes.append(decomp)
                  
        ConnectAtts(posiNodes, "tangentUx", fbfmNodes, "in10")
        ConnectAtts(posiNodes, "tangentUy", fbfmNodes, "in11")
        ConnectAtts(posiNodes, "tangentUz", fbfmNodes, "in12")
        ConnectAtts(posiNodes, "normalX", fbfmNodes, "in00")
        ConnectAtts(posiNodes, "normalY", fbfmNodes, "in01")
        ConnectAtts(posiNodes, "normalZ", fbfmNodes, "in02")
        ConnectAtts(posiNodes, "tangentVx", fbfmNodes, "in20")
        ConnectAtts(posiNodes, "tangentVy", fbfmNodes, "in21")
        ConnectAtts(posiNodes, "tangentVz", fbfmNodes, "in22")
        ConnectAtts(fbfmNodes, "output", decompNodes, "inputMatrix")
        ConnectAtts(decompNodes, "outputRotate", groups, "r")
        ConnectAtts(posiNodes, "position", groups, "t")
        
        if scale is not None:
            ConnectSurfaceToScale(sampleSurface, uvVal, groups, scale, uvScale)


def DoFrameCache(inName="fCache", inNode="locator1", outNodes=["locator9"]):
    # Get objects which trail another object based on time
    data = {}
    mTime = maya.OpenMaya.MTime(1, maya.OpenMaya.MTime.k6000FPS)
    factr = mTime.asUnits(maya.OpenMaya.MTime.uiUnit())
    data["decomp"] = maya.cmds.createNode("decomposeMatrix", n="%s_decomp" % inName, ss=True)
    
    for j in range(len(outNodes)):
        i = j+1
        nodes ={}
        name = "%s_%02d" % (inName, i)
        newAttr = "offset_%02d" % i
        
        nodes["recomp"] = maya.cmds.createNode("composeMatrix", n="%s_recomp" % name, ss=True)
        nodes["compit"] = maya.cmds.createNode("decomposeMatrix", n="%s_compit" % name, ss=True)
        nodes["animBlend"] = maya.cmds.createNode("animBlendNodeTime", n="%s_animBlend" % name, ss=True )
        convr = maya.cmds.createNode("timeToUnitConversion", n="%s_ttuc" % name, ss=True)
        maya.cmds.setAttr("%s.conversionFactor" % convr, factr)
        
        if not maya.cmds.attributeQuery(newAttr, n=inNode, ex=True):
            maya.cmds.addAttr(inNode, at="time", ln=newAttr, k=True)

        maya.cmds.setAttr("%s.%s" % (inNode, newAttr), -i)
        maya.cmds.setAttr("%s.useEulerRotation" % nodes["recomp"],  0)
        maya.cmds.connectAttr("time1.outTime", "%s.inputA" % nodes["animBlend"])
        maya.cmds.connectAttr("%s.%s" % (inNode, newAttr), "%s.inputB" % nodes["animBlend"])
        maya.cmds.connectAttr("%s.output" % nodes["animBlend"], "%s.input" % convr)
        maya.cmds.connectAttr("%s.worldMatrix[0]" % inNode, "%s.inputMatrix" % data["decomp"])
        maya.cmds.connectAttr("%s.outputMatrix" % nodes["recomp"], "%s.inputMatrix" % nodes["compit"])
        
        for node in nodes.values():
            maya.cmds.setAttr("%s.isHistoricallyInteresting" % node, False)
            
        nodes["frameCache"] = {}
        for ts in ["Translate", "Scale"]:
            for xyz in "XYZ":
                nodes["frameCache"]["%s%s" % (ts, xyz)] = maya.cmds.createNode("frameCache",n="%sCache_%s%s" % (name, ts, xyz), ss=True)
        for wxyz in "WXYZ":
            nodes["frameCache"]["Quat%s" % wxyz] = maya.cmds.createNode("frameCache",n="%sCache_r%s" % (name, wxyz), ss=True)
            
        for node in nodes["frameCache"].values():
            maya.cmds.setAttr("%s.isHistoricallyInteresting" % node, False)
        
        for one in nodes["frameCache"]:
            maya.cmds.connectAttr("%s.output" % convr, "%s.varyTime" % nodes["frameCache"][one])
            maya.cmds.connectAttr("%s.output%s" % (data["decomp"], one), "%s.stream" % nodes["frameCache"][one])
            maya.cmds.connectAttr("%s.varying" % nodes["frameCache"][one], "%s.input%s" % (nodes["recomp"], one))
        for trs in ["translate", "rotate", "scale"]:
            maya.cmds.connectAttr("%s.output%s" % (nodes["compit"], trs.capitalize()), "%s.%s" % (outNodes[j], trs))
                
        data[i] = nodes
    return data
    