function [S,GenSys,SVCD,TCSCD,UPFCD,HVDCCD,geib_idx,gsib_idx,lmod_con,rlmod_con,load_con] = rdpstv2_modified(dfile,pathname)
% Reads pst data m file and puts data into class structures.  

% Read PSTV2 data file 
bus = []; 
line = []; exc_con=[];mac_con=[];pss_con = [];hgt_con = [];tgt_con = [];
svc_con = [];tcsc_con = [];ind_con = [];upfc_con = [];
dcl_con = [];dcr_con = [];dci_con = [];
dcld_con = [];dcrc_con = [];dcic_con = [];

wd = cd;
S = InitNetStr;
GenSys=[];HVDCCD=[];geib_idx=[];
gsib_idx=[];lmod_con=[];rlmod_con=[];load_con=[];
options.Resize = 'on';
% [dfile,pathname]=uigetfile('d*.m','Select Data File');
% if pathname == 0
%     error('read pstv2 error, no data file selected')
% else
%     disp(['User selected ', fullfile(pathname, dfile)])
    % cd(pathname);
    lfile =length(dfile);
    % strip off .m and convert to lower case
    %dfile = lower(dfile(1:lfile-2));
    dfile = dfile(1:lfile-2);
    eval(dfile);
    cd(wd);
% end
% check for valid dynamic data file
if isempty(bus)
    error(' the selected file is not a valid data file')
end

% base data
basdat{1} = '100';
basdat{2} = '60';
% basdat = inputdlg({'Base MVA:','Base Frequency Hz:'},'Input Base Data',1,{'100','60'},options); 
%sys_freq = str2num(basdat{2});
%basrad = 2*pi*sys_freq; % default system frequency is 60 Hz
basmva = str2num(basdat{1});
S.BasMva = basmva;


%headers
headers{1} = dfile;
headers{2} = 'base';
% headers = inputdlg({'Header 1:','Header 2:'},'Input Case Headers',1,{dfile,'base'},options);
S.Header1 = headers{1};
S.Header2 = headers{2};

switch size(bus,2)
    case 10
        bus(:,11)= zeros(size(bus,1),1);
        bus(:,12)= zeros(size(bus,1),1);
        bus(:,13)=ones(size(bus,1),1);
        bus(:,14) = 1.5*ones(size(bus,1),1);
        bus(:,15) = 0.5*ones(size(bus,1),1);
    case 12
        bus(:,13)=ones(size(bus,1),1);
        bus(:,14) = 1.5*ones(size(bus,1),1);
        bus(:,15) = 0.5*ones(size(bus,1),1);
    case 13
        bus(:,14) = 1.5*ones(size(bus,1),1);
        bus(:,15) = 0.5*ones(size(bus,1),1);
end
if size(line,2)<8
    line(:,8:10)=zeros(size(line,1),3);
end

% input bus data
n_bus = size(bus,1);
Bus.NumBus = n_bus; Bus.BaseVolts = bus(:,13);
Bus.BusNum = cellstr(int2str(bus(:,1)));
Bus.BusName = cellstr(int2str(bus(:,1)));
Bus.Type = bus(:,10);
Bus.G = S.BasMva*sparse(bus(:,8)); Bus.B = S.BasMva*sparse(bus(:,9));
Bus.Area=ones(Bus.NumBus,1);Bus.Zone = Bus.Area;
Bus.VoltMag = bus(:,2); Bus.VoltAngd = bus(:,3);
Bus.VoltAngr = Bus.VoltAngd.*pi/180; 
Bus.Owner = Bus.Area;
bn = bus(:,1);
Bus.Index = sparse(max(bn),1); Bus.Index(bn)=(1:n_bus)';
S.Bus = Bus;
% input load data
lb_idx = find(bus(:,6)|bus(:,7));
NumLoads = length(lb_idx);
Status = ones(NumLoads,1);
Load.BusRef = cellstr(num2str(bus(lb_idx,1)));
Load.MultLoadId = cellstr(int2str(ones(n_bus,1))); 
Load.Status = Status;Load.Area = Status;Load.Zone = Status;
Load.ConstP = S.BasMva*bus(lb_idx,6);
Load.ConstQ = S.BasMva*bus(lb_idx,7);
Load.ConstIP = sparse(NumLoads,1); Load.ConstIQ = sparse(NumLoads,1);
Load.ConstZP = sparse(NumLoads,1);Load.ConstZQ = sparse(NumLoads,1);
Load.Owner = Status;
Load.Index = lb_idx;
S.Load = Load; clear Load;

% input Generator data
gb_idx = find(bus(:,10)<3);
n_gen = length(gb_idx);
Gen.BusRef = cellstr(int2str(bus(gb_idx,1))) ;
Gen.MultGenId = cellstr(num2str(ones(n_gen,1)));
Gen.PGen = S.BasMva*bus(gb_idx,4);Gen.QGen = S.BasMva*bus(gb_idx,5);
Gen.QMax = S.BasMva*bus(gb_idx,11); Gen.QMin = S.BasMva*bus(gb_idx,12);
zqm = find(Gen.QMax==0&Gen.QMin==0);
if~isempty(zqm)
    Gen.QMax(zqm)=S.BasMva*999*ones(length(zqm),1);
    Gen.QMin(zqm)=-S.BasMva*999*ones(length(zqm),1);
end
Gen.PMax = 1e3*S.BasMva*ones(n_gen,1); Gen.PMin = -1e6*0*S.BasMva*ones(n_gen,1); Gen.VSet = bus(gb_idx,2);
Gen.RegBus = Gen.BusRef; Gen.RSource = zeros(n_gen,1); Gen.XSource =zeros(n_gen,1);
Gen.RTF =zeros(n_gen,1); Gen.XTF = zeros(n_gen,1); Gen.Tap = ones(n_gen,1); 
Gen.Status = ones(n_gen,1);Gen.ContProp = 100*ones(n_gen,1);
Gen.BasMva = S.BasMva*ones(n_gen,1);
Gen.Owner = []; Gen.OwnerProp = [];
Gen.BusIndex = gb_idx;
Gen.RegBusIndex = gb_idx;Gen.RegIndex = [];
S.Gen = Gen; clear Gen;

% input the branch data
n_line = length(line(:,1));
Line.FromBusRef = cellstr(int2str(line(:,1))); Line.ToBusRef = cellstr(int2str(line(:,2)));
fi = com_index(bus(:,1),line(:,1));
Line.FromIndex = fi;
ti = com_index(bus(:,1),line(:,2));
Line.ToIndex = ti;
for k = 1:n_line
    ml_idx = find(line(:,1)==line(k,1)& line(:,2)==line(k,2));
    for kk = 1:length(ml_idx)
        Line.MultCctId(ml_idx(kk),1)={num2str(kk)};
    end
end
Line.R = line(:,3); Line.X = line(:,4); Line.Bc = line(:,5);
Line.CurRat1 = zeros(n_line,1);Line.CurRat2 = zeros(n_line,1);Line.CurRat3 = zeros(n_line,1);
Line.Tap = line(:,6); Line.Phase = line(:,7); 
Line.FromG = zeros(n_line,1);Line.FromB=zeros(n_line,1);
Line.ToG = zeros(n_line,1);Line.ToB=zeros(n_line,1);
Line.Status = ones(n_line,1); Line.Length = zeros(n_line,1); 
Line.Owner = ones(n_line,4); 
Line.OwnerProp = ones(n_line,4);
Line.Meter = sparse(n_line,1);

S.Line = Line;
% Get the transformer control data 
nlc = size(line,2);
if nlc<10
    % no transformer
    Tran = [];
else
    tran_idx = find(line(:,8)~=0);
    if ~isempty(tran_idx)
        n_tran = length(tran_idx);
        Tran.LineIndex = tran_idx;
        for k = 1:n_tran
            Tran.FromBusRef(k,1) = {char(Line.FromBusRef(tran_idx(k)))}; 
            Tran.ToBusRef(k,1) = {char(Line.ToBusRef(tran_idx(k)))}; 
        end
        Tran.FromIndex = Line.FromIndex(tran_idx);
        Tran.ToIndex = Line.ToIndex(tran_idx);
        Tran.Meter = zeros(n_tran,1);% controlled bus assummed to be at to bus side
        Tran.MultCctId = Line.MultCctId(tran_idx); 
        Tran.ContRef = Tran.ToBusRef;
        Tran.ContIndex=Tran.ToIndex;
        Tran.MaxTap = line(tran_idx,8);
        Tran.MinTap = line(tran_idx,9);
        Tran.MaxV = bus(Line.ToIndex(tran_idx),14);
        Tran.MinV = bus(Line.ToIndex(tran_idx),15);
        Tran.Step = line(tran_idx,10); 
        Tran.TabNum = zeros(n_tran,1); 
        Tran.Freeze = zeros(n_tran,1);
        Tran.RComp = zeros(n_tran,1);
        Tran.XComp = zeros(n_tran,1);
    else
        Tran = [];
    end
end
S.Tran = Tran; clear Tran;clear Line line
% input induction motor data
if ~isempty(ind_con)
	S.IND.MacNum = ind_con(:,1);
	S.IND.BusNum = ind_con(:,2);
    MotBusIndex = full(S.Bus.Index(ind_con(:,2)));
    S.IND.BusIndex = MotBusIndex;
	S.IND.BasMva = ind_con(:,3);
	S.IND.rs = ind_con(:,4);
	S.IND.xs = ind_con(:,5);
	S.IND.Xm = ind_con(:,6);
	S.IND.rr  = ind_con(:,7);
	S.IND.xr = ind_con(:,8);
	S.IND.rr1 = ind_con(:,10);
	S.IND.xr1 = ind_con(:,11);
	S.IND.H = ind_con(:,9);
	S.IND.dbf = ind_con(:,12);
	S.IND.ISat = ind_con(:,13);
    S.IND.kst = mld_con(ind_con(:,1),3);
    S.IND.tsti = mld_con(ind_con(:,1),4);
    S.IND.kspd = mld_con(ind_con(:,1),5);
    S.IND.tspdi = mld_con(ind_con(:,1),6);
	LoadBusIndex =S.Load.Index ;
	ml_idx = com_index(LoadBusIndex,MotBusIndex);
	S.IND.PElec = S.Load.ConstP(ml_idx).*ind_con(:,15);% zero ind_con(:,15) indicates starting motor
	S.Load.ConstP(ml_idx)=S.Load.ConstP(ml_idx)-S.IND.PElec;% reduce power load by machine power
    S.IND.TermVolts = S.Bus.VoltMag(MotBusIndex).*exp(i*S.Bus.VoltAngr(MotBusIndex));
else
	S.IND =[];
end
S.Areas = [];

if isempty(dcl_con)||isempty(dcr_con)||isempty(dci_con)
	S.DC2 = [];
else
	n_dc = length(dcl_con(:,1));
	S.DC2.Line.LineNum = dcl_con(:,1);
	S.DC2.Line.ContMode = dcl_con(:,3);
	S.DC2.Line.LineRes = dcl_con(:,2);
	S.DC2.Line.ContVal = dcl_con(:,4);
	S.DC2.Line.SchedVolts = dcl_con(:,5);
	S.DC2.Line.SwitchVolts = dcl_con(:,6);
	S.DC2.Line.RComp = dcl_con(:,7);
	S.DC2.Line.CurMarg = dcl_con(:,8);
	S.DC2.Line.Meter = [];
	S.DC2.Line.VoltMin = dcl_con(:,9);
    % col 1		Number of Bridges
    % col 2 	alphamax (rectifier) gammamax (inverter)(degrees)
    % col 3 	alphamin (rectifier) gammamin (inverter) (degrees)
    % col 4 	XCom commutating reactance per bridge (ohm)
    % col 5 	HT Base Voltage (KV)
    % col 6     LT Base Voltage (KV) 
    % col 7 	Tap position
    % col 8 	TapMax
    % col 9 	TapMin
    % col 10 	TapStep
    % col 11 	HTBusNumber
    % col 12	LTBusNumber
	S.DC2.Rectifier.BusRef = cellstr(num2str(1:length(dcr_con(:,1))));
	S.DC2.Rectifier.NumBridge = dcr_con(:,1);
	S.DC2.Rectifier.AngMax = dcr_con(:,2);
	S.DC2.Rectifier.AngMin = dcr_con(:,3);
	S.DC2.Rectifier.RCom = zeros(n_dc,1);
	S.DC2.Rectifier.XCom = dcr_con(:,4);
	S.DC2.Rectifier.BaseVolts = dcr_con(:,5);
    S.DC2.Rectifier.TTRatio = dcr_con(:,7);
	S.DC2.Rectifier.Tap = ones(n_dc,1);
	S.DC2.Rectifier.TapMax = dcr_con(:,8);
	S.DC2.Rectifier.TapMin = dcr_con(:,9);
	S.DC2.Rectifier.TapStep = dcr_con(:,10);
	S.DC2.Rectifier.AlphaMeaBus = zeros(n_dc,1);
	S.DC2.Rectifier.FromBus = zeros(n_dc,1);
	S.DC2.Rectifier.ToBus = zeros(n_dc,1);
	S.DC2.Rectifier.CctId = zeros(n_dc,1);
	S.DC2.Rectifier.XCap = zeros(n_dc,1); 
	S.DC2.Rectifier.HTBusIndex = full(Bus.Index(dcr_con(:,11)));
	S.DC2.Rectifier.HTBusName = num2str(dcr_con(:,11));
    ltb = dcr_con(:,12);
    ltnz = find(ltb);
    S.DC2.Rectifier.LTBusName = num2str(ltb);
    S.DC2.Rectifier.LTBusIndex = zeros(n_dc,1);
    if ~isempty(ltnz)
        S.DC2.Rectifier.LTBusIndex(ltnz) = full(Bus.Index(dci_con(ltnz,12)));
    end
	S.DC2.Inverter.BusRef = cellstr(num2str(1:length(dci_con(:,1))));
	S.DC2.Inverter.NumBridge = dci_con(:,1);
	S.DC2.Inverter.AngMax = dci_con(:,2);
	S.DC2.Inverter.AngMin = dci_con(:,3);
	S.DC2.Inverter.RCom = zeros(n_dc,1);
	S.DC2.Inverter.XCom = dci_con(:,4);
	S.DC2.Inverter.BaseVolts = dci_con(:,5);
	S.DC2.Inverter.TTRatio = dci_con(:,7);
    S.DC2.Inverter.Tap = ones(n_dc,1);
	S.DC2.Inverter.TapMax = dci_con(:,8);
	S.DC2.Inverter.TapMin = dci_con(:,9);
	S.DC2.Inverter.TapStep = dci_con(:,10);
	S.DC2.Inverter.AlphaMeaBus = zeros(n_dc,1);
	S.DC2.Inverter.FromBus = zeros(n_dc,1);
	S.DC2.Inverter.ToBus = zeros(n_dc,1);
	S.DC2.Inverter.CctId = zeros(n_dc,1);
	S.DC2.Inverter.XCap = zeros(n_dc,1);
	S.DC2.Inverter.HTBusIndex = full(Bus.Index(dci_con(:,11)));
	S.DC2.Inverter.HTBusName = num2str(dci_con(:,11));
	ltb = dci_con(:,12);
    ltnz = find(ltb);
    S.DC2.Inverter.LTBusName = num2str(ltb);
    S.DC2.Inverter.LTBusIndex = zeros(n_dc,1);
    if ~isempty(ltnz)
        S.DC2.Inverter.LTBusIndex(ltnz) = full(Bus.Index(dci_con(ltnz,12)));
    end
end

if ~isempty(svc_con)
    sbn = svc_con(:,2);
    n_svc = length(sbn);
    sbi = com_index(bus(:,1),svc_con(:,2));
    Bs = svc_con(:,10).*bus(sbi,9);% shunt B on system base
    S.Shunt.BusNum = svc_con(:,2);
    S.Shunt.BusIndex = sbi;
    S.Shunt.BusRef = cellstr(num2str(svc_con(:,2)));
    S.Shunt.ContBusIndex = sbi;
    S.Shunt.Mode = 2*ones(n_svc,1);
    S.Shunt.VMax = S.Bus.VoltMag(sbi);
    S.Shunt.VMin = S.Bus.VoltMag(sbi);
    S.Shunt.ContBusRef=S.Shunt.BusRef;
    S.Bus.B(sbi)=S.Bus.B(sbi)-S.BasMva*Bs;
    S.Shunt.BStart = S.BasMva*Bs;
    S.Shunt.NumStep = sparse(n_svc,8);
    S.Shunt.BVal = sparse(n_svc,8);
    S.Shunt.NumStep(:,1) = 1;
    S.Shunt.BVal(:,1) = svc_con(:,5).*svc_con(:,3);
    S.Shunt.NumStep(:,2) = 1;
    S.Shunt.BVal(:,2) = svc_con(:,4).*svc_con(:,3);
else
    S.Shunt = [];
end
if ~isempty(upfc_con)
	upfcn = upfc_con(:,1);
	n_upfc =length(upfcn);
    S.FACTS.Num = upfc_con(:,1);
	S.FACTS.SendBus = upfc_con(:,2);
	S.FACTS.TermEndBus = upfc_con(:,3);
	S.FACTS.Mode = upfc_con(:,4);
	S.FACTS.P = basmva*upfc_con(:,5);
	S.FACTS.Q = basmva*upfc_con(:,6);
	S.FACTS.VSet =upfc_con(:,7);
	S.FACTS.SCMax = basmva*upfc_con(:,8);
	S.FACTS.BTMax = basmva*upfc_con(:,9);
	S.FACTS.VTMin = upfc_con(:,10);
	S.FACTS.VTMax =upfc_con(:,11);
	S.FACTS.CSMax = basmva*upfc_con(:,12);
	S.FACTS.Xs = upfc_con(:,13);
	S.FACTS.Xsh = upfc_con(:,14);
	S.FACTS.Owner = ones(n_upfc,1);
else
	S.FACTS = [];
end
clear bus
S.TTab = [];

S.MTDC = [];

S.MSLG = [];

S.Zone = [];

S.Atd = [];


S.Owner = [];

GenSys.Gen = [];GenSys.EXCDC = [];GenSys.EXCAC = []; GenSys.EXCST = [];GenSys.HGovD = [];GenSys.TGovD = []; 
SVCD = [];TCSCD = [];UPFCD=[];
if ~isempty(mac_con)
    n_gen = length(mac_con(:,1));
    GenSys.Gen.NumGen = n_gen;
    GenSys.Gen.Type = mac_con(:,19);
    GenSys.Gen.Num = mac_con(:,1);
    GenSys.Gen.Name = cellstr(num2str(mac_con(:,1)));
    GenSys.Gen.BusNum = mac_con(:,2);
    GenSys.Gen.BasMva = mac_con(:,3);
    GenSys.Gen.xl = mac_con(:,4);
    GenSys.Gen.ra = mac_con(:,5);
    GenSys.Gen.xd = mac_con(:,6);
    GenSys.Gen.xdp = mac_con(:,7);
    GenSys.Gen.xdpp = mac_con(:,8);
    GenSys.Gen.Tdpo = mac_con(:,9);
    GenSys.Gen.Tdppo = mac_con(:,10);
    GenSys.Gen.xq = mac_con(:,11);
    GenSys.Gen.xqp =mac_con(:,12);
    GenSys.Gen.xqpp = mac_con(:,13);
    GenSys.Gen.Tqpo = mac_con(:,14);
    GenSys.Gen.Tqppo = mac_con(:,15);
    GenSys.Gen.H = mac_con(:,16);
    GenSys.Gen.damp = mac_con(:,17);
    GenSys.Gen.S1 = mac_con(:,20);
    GenSys.Gen.S2 = mac_con(:,21);
    if size(mac_con,2)==23
        mc22z = find(mac_con(:,22)==0);
        mc23z = find(mac_con(:,23)==0);
        GenSys.Gen.apf = mac_con(:,22);
        GenSys.Gen.rpf = mac_con(:,23);
        GenSys.Gen.apf(mc22z) = 1;
        GenSys.Gen.rpf(mc23z) = 1;
    else
        GenSys.Gen.apf = ones(n_gen,1);
        GenSys.Gen.rpf = ones(n_gen,1);
    end
    if ~isempty(exc_con)
        n_exc = length(exc_con(:,1));
        [bexc,iexc] = unique(exc_con(:,2));
        if length(bexc)~=n_exc
            uiwait(msgbox('there are multiple exciters at some generators','exciter data error','modal'))
            exc_con=exc_con(iexc,:);
            %n_exc = length(exc_con(:,2));
        end
        Type = exc_con(:,1);
        dc12_idx = find(Type==1|Type==2);
        dc3_idx = find(Type == 3);
        %dc_idx = sort([dc12_idx;dc3_idx]);
        %dc12_idx = com_index(dc_idx,dc12_idx);
        %dc3_idx = com_index(dc_idx,dc3_idx);
        ac1_idx = find(Type==4);
        ac2_idx = find(Type==5);
        ac3_idx = find(Type==6);
        ac4_idx = find(Type==7);
        ac5_idx = find(Type==8);
        ac6_idx = find(Type==9);
        ac7_idx = find(Type==13);
        ac8_idx = find(Type==16);
        
        st1_idx = find(Type==10);
        st2_idx = find(Type==11);
        st3_idx = find(Type==12);
        st4_idx = find(Type==14);
        st5_idx = find(Type==15);
        

        Rc = exc_con(:,3);
        Xc = exc_con(:,4); 
        Tr = exc_con(:,5); 
        if ~isempty(pss_con)
            PGenNum = pss_con(:,2);
        else
            PGenNum = [];
        end
        GenSys.EXCDC.e12 = [];GenSys.EXCDC.a12 = [];GenSys.EXCDC.e3 = [];GenSys.EXCDC.a3 = [];
        if ~isempty(dc12_idx)
            GenNum = exc_con(dc12_idx,2); 
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                PSSD = pdata(pss_con,pss_idx);
                GenSys.EXCDC.a12.PSSD = PSSD;
            else 
                GenSys.EXCDC.a12.PSSD = [];
            end
            GenSys.EXCDC.a12.TransD.Type = Type(dc12_idx);
            GenSys.EXCDC.a12.TransD.GenNum = GenNum;
            GenSys.EXCDC.a12.TransD.ExcIndex = dc12_idx;
            GenSys.EXCDC.a12.TransD.Rc = Rc(dc12_idx);
            GenSys.EXCDC.a12.TransD.Xc = Xc(dc12_idx);
            GenSys.EXCDC.a12.TransD.Tr = Tr(dc12_idx);
            GenSys.EXCDC.e12.Type = Type(dc12_idx);
            GenSys.EXCDC.e12.GenNum = GenNum;
            GenSys.EXCDC.e12.ExcIndex = dc12_idx;
            GenSys.EXCDC.e12.VrMx = exc_con(dc12_idx,18);
%             GenSys.EXCDC.e12.VrMn = exc_con(dc12_idx,19);
%             GenSys.EXCDC.e12.Ke = exc_con(dc12_idx,21);
%             GenSys.EXCDC.e12.Te = exc_con(dc12_idx,22);
%             GenSys.EXCDC.e12.se1 = exc_con(dc12_idx,23);
%             GenSys.EXCDC.e12.se2 = exc_con(dc12_idx,24);
%             GenSys.EXCDC.e12.ve1 = exc_con(dc12_idx,25);
%             GenSys.EXCDC.e12.ve2 = exc_con(dc12_idx,26);
            GenSys.EXCDC.a12.GenNum = exc_con(dc12_idx,2);
            GenSys.EXCDC.a12.ExcIndex = dc12_idx;
            GenSys.EXCDC.a12.Ka = exc_con(dc12_idx,6);
            GenSys.EXCDC.a12.Ta = exc_con(dc12_idx,7);
            GenSys.EXCDC.a12.Tb = exc_con(dc12_idx,8);
            GenSys.EXCDC.a12.Tc = exc_con(dc12_idx,9);
            GenSys.EXCDC.a12.Kf = exc_con(dc12_idx,10);
            GenSys.EXCDC.a12.Tf = exc_con(dc12_idx,11);
        end
        if ~isempty(dc3_idx)
            % pss not possible on dc3 exciters
            GenSys.EXCDC.a3.TransD.Type = Type(dc3_idx);
            GenSys.EXCDC.a3.TransD.GenNum = exc_con(dc3_idx,2);
            GenSys.EXCDC.a3.TransD.ExcIndex = dc3_idx; 
            GenSys.EXCDC.a3.TransD.Rc = Rc(dc3_idx);
            GenSys.EXCDC.a3.TransD.Xc = Xc(dc3_idx);
            GenSys.EXCDC.a3.TransD.Tr = Tr(dc3_idx);
            GenSys.EXCDC.e3.Type = Type(dc3_idx);
            GenSys.EXCDC.e3.GenNum = exc_con(dc3_idx,2);
            GenSys.EXCDC.e3.ExcIndex = dc3_idx;
            GenSys.EXCDC.e3.VrMx = exc_con(dc3_idx,18);
            GenSys.EXCDC.e3.VrMn = exc_con(dc3_idx,19);
            GenSys.EXCDC.e3.Ke = exc_con(dc3_idx,21);
            GenSys.EXCDC.e3.Te = exc_con(dc3_idx,22);
            GenSys.EXCDC.e3.se1 = exc_con(dc3_idx,23);
            GenSys.EXCDC.e3.se2 = exc_con(dc3_idx,24);
            GenSys.EXCDC.e3.ve1 = exc_con(dc3_idx,25);
            GenSys.EXCDC.e3.ve2 = exc_con(dc3_idx,26);
            GenSys.EXCDC.a3.GenNum = exc_con(dc3_idx,2);
            GenSys.EXCDC.a3.ExcIndex = dc3_idx;
            Kv = exc_con(dc3_idx,6);
            Kvz_idx = find(Kv==0);
            if ~isempty(Kvz_idx)
                Kv(Kvz_idx) = 1;
                uiwait(msgbox('Zero Kv set to 1','dc3a data error','modal'))
            end
            Trh = exc_con(dc3_idx,7);
            Trhz_idx = find(Trh==0);
            if ~isempty(Trhz_idx)
                Trh(Trhz_idx) = 20;
                uiwait(msgbox('Zero Trh set to 20','dc3a data error','modal'))
            end
            GenSys.EXCDC.a3.Kv = Kv;
            GenSys.EXCDC.a3.Trh = Trh;
        end
        GenSys.EXCAC.e1 = [];GenSys.EXCAC.a1 = [];
        GenSys.EXCAC.e2 = [];GenSys.EXCAC.a2 = [];
        GenSys.EXCAC.e3 = [];GenSys.EXCAC.a3 = [];
        GenSys.EXCAC.e4 = [];GenSys.EXCAC.a4 = [];
        GenSys.EXCAC.e5 = [];GenSys.EXCAC.a5 = [];
        GenSys.EXCAC.e6 = [];GenSys.EXCAC.a6 = [];
        GenSys.EXCAC.e7 = [];GenSys.EXCAC.a7 = [];
        GenSys.EXCAC.e8 = [];GenSys.EXCAC.a8 = [];
        if ~isempty(ac1_idx)
            GenNum = exc_con(ac1_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a1.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a1.PSSD = [];
            end
            GenSys.EXCAC.e1.Type = Type(ac1_idx);
            GenSys.EXCAC.e1.GenNum = GenNum;
            GenSys.EXCAC.e1.ExcIndex = ac1_idx;
            GenSys.EXCAC.e1.Kc = exc_con(ac1_idx,20);
            GenSys.EXCAC.e1.Ke = exc_con(ac1_idx,21);
            GenSys.EXCAC.e1.Te = exc_con(ac1_idx,22);
            GenSys.EXCAC.e1.se1 = exc_con(ac1_idx,23);
            GenSys.EXCAC.e1.se2 = exc_con(ac1_idx,24);
            GenSys.EXCAC.e1.ve1 = exc_con(ac1_idx,25);
            GenSys.EXCAC.e1.ve2 = exc_con(ac1_idx,26);
            GenSys.EXCAC.e1.Kd = exc_con(ac1_idx,29);
            GenSys.EXCAC.a1.TransD.Type = Type(ac1_idx);
            GenSys.EXCAC.a1.TransD.GenNum = exc_con(ac1_idx,2);
            GenSys.EXCAC.a1.TransD.ExcIndex = ac1_idx; 
            GenSys.EXCAC.a1.TransD.Rc = Rc(ac1_idx);
            GenSys.EXCAC.a1.TransD.Xc = Xc(ac1_idx);
            GenSys.EXCAC.a1.TransD.Tr = Tr(ac1_idx);
            GenSys.EXCAC.a1.Type = Type(ac1_idx);
            GenSys.EXCAC.a1.GenNum = exc_con(ac1_idx,2);
            GenSys.EXCAC.a1.ExcIndex = ac1_idx;
            GenSys.EXCAC.a1.Ka = exc_con(ac1_idx,6);
            GenSys.EXCAC.a1.Ta = exc_con(ac1_idx,7);
            GenSys.EXCAC.a1.Tb = exc_con(ac1_idx,8);
            GenSys.EXCAC.a1.Tc = exc_con(ac1_idx,9);
            GenSys.EXCAC.a1.Kf = exc_con(ac1_idx,10);
            GenSys.EXCAC.a1.Tf = exc_con(ac1_idx,11);
            GenSys.EXCAC.a1.VrMx = exc_con(ac1_idx,18);
            GenSys.EXCAC.a1.VrMn = exc_con(ac1_idx,19);
            GenSys.EXCAC.a1.VaMx = exc_con(ac1_idx,16);
            GenSys.EXCAC.a1.VaMn = exc_con(ac1_idx,17);
        end
        if~isempty(ac2_idx)
            GenNum = exc_con(ac2_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a2.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a2.PSSD = [];
            end
            GenSys.EXCAC.e2.Type = Type(ac2_idx);
            GenSys.EXCAC.e2.GenNum = GenNum;
            GenSys.EXCAC.e2.ExcIndex = ac2_idx;
            GenSys.EXCAC.e2.Kc = exc_con(ac2_idx,20);
            GenSys.EXCAC.e2.Ke = exc_con(ac2_idx,21);
            GenSys.EXCAC.e2.Te = exc_con(ac2_idx,22);
            GenSys.EXCAC.e2.se1 = exc_con(ac2_idx,23);
            GenSys.EXCAC.e2.se2 = exc_con(ac2_idx,24);
            GenSys.EXCAC.e2.ve1 = exc_con(ac2_idx,25);
            GenSys.EXCAC.e2.ve2 = exc_con(ac2_idx,26);
            GenSys.EXCAC.e2.Kd = exc_con(ac2_idx,29);
            GenSys.EXCAC.e2.Vfemx = exc_con(ac2_idx,30);
            GenSys.EXCAC.a2.TransD.Type = Type(ac2_idx);
            GenSys.EXCAC.a2.TransD.GenNum = exc_con(ac2_idx,2);
            GenSys.EXCAC.a2.TransD.ExcIndex = ac2_idx; 
            GenSys.EXCAC.a2.TransD.Rc = Rc(ac2_idx);
            GenSys.EXCAC.a2.TransD.Xc = Xc(ac2_idx);
            GenSys.EXCAC.a2.TransD.Tr = Tr(ac2_idx);
            GenSys.EXCAC.a2.Type = Type(ac2_idx);
            GenSys.EXCAC.a2.GenNum = exc_con(ac2_idx,2);
            GenSys.EXCAC.a2.ExcIndex = ac2_idx;
            GenSys.EXCAC.a2.Ka = exc_con(ac2_idx,6);
            GenSys.EXCAC.a2.Ta = exc_con(ac2_idx,7);
            GenSys.EXCAC.a2.Tb = exc_con(ac2_idx,8);
            GenSys.EXCAC.a2.Tc = exc_con(ac2_idx,9);
            GenSys.EXCAC.a2.Kf = exc_con(ac2_idx,10);
            GenSys.EXCAC.a2.Tf = exc_con(ac2_idx,11);
            GenSys.EXCAC.a2.VrMx = exc_con(ac2_idx,18);
            GenSys.EXCAC.a2.VrMn = exc_con(ac2_idx,19);
            GenSys.EXCAC.a2.VaMx = exc_con(ac2_idx,16);
            GenSys.EXCAC.a2.VaMn = exc_con(ac2_idx,17);
            GenSys.EXCAC.a2.Kh = exc_con(ac2_idx,12);
            GenSys.EXCAC.a2.Kb = exc_con(ac2_idx,13);
        end
        if~isempty(ac3_idx)
            n_ac3 = length(ac3_idx);
            GenNum = exc_con(ac3_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a3.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a3.PSSD = [];
            end
            GenSys.EXCAC.e3.Type = Type(ac3_idx);
            GenSys.EXCAC.e3.GenNum = GenNum;
            GenSys.EXCAC.e3.ExcIndex = ac3_idx;
            GenSys.EXCAC.e3.Kc = exc_con(ac3_idx,20);
            GenSys.EXCAC.e3.Ke = exc_con(ac3_idx,21);
            GenSys.EXCAC.e3.Te = exc_con(ac3_idx,22);
            GenSys.EXCAC.e3.se1 = exc_con(ac3_idx,23);
            GenSys.EXCAC.e3.se2 = exc_con(ac3_idx,24);
            GenSys.EXCAC.e3.ve1 = exc_con(ac3_idx,25);
            GenSys.EXCAC.e3.ve2 = exc_con(ac3_idx,26);
            GenSys.EXCAC.e3.Kd = exc_con(ac3_idx,29);
            GenSys.EXCAC.e3.Vfemx = exc_con(ac3_idx,30);
            GenSys.EXCAC.e3.VeMn = sparse(n_ac3,1);
            GenSys.EXCAC.e3.Kr = exc_con(ac3_idx,27);
            GenSys.EXCAC.a3.TransD.Type = Type(ac3_idx);
            GenSys.EXCAC.a3.TransD.GenNum = exc_con(ac3_idx,2);
            GenSys.EXCAC.a3.TransD.ExcIndex = ac3_idx; 
            GenSys.EXCAC.a3.TransD.Rc = Rc(ac3_idx);
            GenSys.EXCAC.a3.TransD.Xc = Xc(ac3_idx);
            GenSys.EXCAC.a3.TransD.Tr = Tr(ac3_idx);
            GenSys.EXCAC.a3.Type = Type(ac3_idx);
            GenSys.EXCAC.a3.GenNum = exc_con(ac3_idx,2);
            GenSys.EXCAC.a3.ExcIndex = ac3_idx;
            GenSys.EXCAC.a3.Ka = exc_con(ac3_idx,6);
            GenSys.EXCAC.a3.Ta = exc_con(ac3_idx,7);
            GenSys.EXCAC.a3.Tb = exc_con(ac3_idx,8);
            GenSys.EXCAC.a3.Tc = exc_con(ac3_idx,9);
            GenSys.EXCAC.a3.Kf = exc_con(ac3_idx,10);
            GenSys.EXCAC.a3.Tf = exc_con(ac3_idx,11);
            GenSys.EXCAC.a3.Efdn = exc_con(ac3_idx,12);
            GenSys.EXCAC.a3.Kn = exc_con(ac3_idx,13);
            GenSys.EXCAC.a3.VrMx = exc_con(ac3_idx,18);
            GenSys.EXCAC.a3.VrMn = exc_con(ac3_idx,19);
        end	
        if~isempty(ac4_idx)
            GenNum = exc_con(ac4_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a4.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a4.PSSD = [];
            end
            GenSys.EXCAC.e4.Type = Type(ac4_idx);
            GenSys.EXCAC.e4.GenNum = GenNum;
            GenSys.EXCAC.e4.ExcIndex = ac4_idx;
            GenSys.EXCAC.e4.Kc = exc_con(ac4_idx,20);
            GenSys.EXCAC.a4.TransD.Type = Type(ac4_idx);
            GenSys.EXCAC.a4.TransD.GenNum = exc_con(ac4_idx,2);
            GenSys.EXCAC.a4.TransD.ExcIndex = ac4_idx; 
            GenSys.EXCAC.a4.TransD.Rc = Rc(ac4_idx);
            GenSys.EXCAC.a4.TransD.Xc = Xc(ac4_idx);
            GenSys.EXCAC.a4.TransD.Tr = Tr(ac4_idx);
            GenSys.EXCAC.a4.Type = Type(ac4_idx);
            GenSys.EXCAC.a4.GenNum = exc_con(ac4_idx,2);
            GenSys.EXCAC.a4.ExcIndex = ac4_idx;
            GenSys.EXCAC.a4.Ka = exc_con(ac4_idx,6);
            GenSys.EXCAC.a4.Ta = exc_con(ac4_idx,7);
            GenSys.EXCAC.a4.Tb = exc_con(ac4_idx,8);
            GenSys.EXCAC.a4.Tc = exc_con(ac4_idx,9);
            GenSys.EXCAC.a4.ViMx = exc_con(ac4_idx,14);
            GenSys.EXCAC.a4.ViMn = exc_con(ac4_idx,15);
            GenSys.EXCAC.a4.VrMx = exc_con(ac4_idx,18);
            GenSys.EXCAC.a4.VrMn = exc_con(ac4_idx,19);
        end
        if~isempty(ac5_idx)
            n_ac5 = length(ac5_idx);
            GenNum = exc_con(ac5_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a5.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a5.PSSD = [];
            end
            GenSys.EXCAC.e5.Type = Type(ac5_idx);
            GenSys.EXCAC.e5.GenNum = GenNum;
            GenSys.EXCAC.e5.ExcIndex = ac5_idx;
            GenSys.EXCAC.e5.Kc = exc_con(ac5_idx,20);
            GenSys.EXCAC.e5.Ke = exc_con(ac5_idx,21);
            GenSys.EXCAC.e5.Te = exc_con(ac5_idx,22);
            GenSys.EXCAC.e5.se1 = exc_con(ac5_idx,23);
            GenSys.EXCAC.e5.se2 = exc_con(ac5_idx,24);
            GenSys.EXCAC.e5.ve1 = exc_con(ac5_idx,25);
            GenSys.EXCAC.e5.ve2 = exc_con(ac5_idx,26);
            GenSys.EXCAC.e5.VeMn = zeros(n_ac5,1);
            GenSys.EXCAC.a5.TransD.Type = Type(ac5_idx);
            GenSys.EXCAC.a5.TransD.GenNum = exc_con(ac5_idx,2);
            GenSys.EXCAC.a5.TransD.ExcIndex = ac5_idx; 
            GenSys.EXCAC.a5.TransD.Rc = Rc(ac5_idx);
            GenSys.EXCAC.a5.TransD.Xc = Xc(ac5_idx);
            GenSys.EXCAC.a5.TransD.Tr = Tr(ac5_idx);
            GenSys.EXCAC.a5.Type = Type(ac5_idx);
            GenSys.EXCAC.a5.GenNum = exc_con(ac5_idx,2);
            GenSys.EXCAC.a5.ExcIndex = ac5_idx;
            GenSys.EXCAC.a5.Ka = exc_con(ac5_idx,6);
            GenSys.EXCAC.a5.Ta = exc_con(ac5_idx,7);
            GenSys.EXCAC.a5.Kf = exc_con(ac5_idx,10);
            GenSys.EXCAC.a5.Tf = exc_con(ac5_idx,11);
            GenSys.EXCAC.a5.Tf1 = exc_con(ac5_idx,12);
            GenSys.EXCAC.a5.Tf2 = exc_con(ac5_idx,13);
            GenSys.EXCAC.a5.VrMx = exc_con(ac5_idx,18);
            GenSys.EXCAC.a5.VrMn = exc_con(ac5_idx,19);
        end	
        if~isempty(ac6_idx)
            n_ac6 = length(ac6_idx);
            GenNum = exc_con(ac6_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a6.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a6.PSSD = [];
            end
            GenSys.EXCAC.e6.Type = Type(ac6_idx);
            GenSys.EXCAC.e6.GenNum = exc_con(ac6_idx,2);
            GenSys.EXCAC.e6.ExcIndex = ac6_idx;
            GenSys.EXCAC.e6.Kc = exc_con(ac6_idx,20);
            GenSys.EXCAC.e6.Ke = exc_con(ac6_idx,21);
            GenSys.EXCAC.e6.Te = exc_con(ac6_idx,22);
            GenSys.EXCAC.e6.se1 = exc_con(ac6_idx,23);
            GenSys.EXCAC.e6.se2 = exc_con(ac6_idx,24);
            GenSys.EXCAC.e6.ve1 = exc_con(ac6_idx,25);
            GenSys.EXCAC.e6.ve2 = exc_con(ac6_idx,26);
            GenSys.EXCAC.e6.Kd = exc_con(ac6_idx,29);
            GenSys.EXCAC.e6.VeMn = zeros(n_ac6,1);
            GenSys.EXCAC.a6.Vfemx = exc_con(ac6_idx,30);
            GenSys.EXCAC.a6.TransD.Type = Type(ac6_idx);
            GenSys.EXCAC.a6.TransD.GenNum = exc_con(ac6_idx,2);
            GenSys.EXCAC.a6.TransD.ExcIndex = ac6_idx; 
            GenSys.EXCAC.a6.TransD.Rc = Rc(ac6_idx);
            GenSys.EXCAC.a6.TransD.Xc = Xc(ac6_idx);
            GenSys.EXCAC.a6.TransD.Tr = Tr(ac6_idx);
            GenSys.EXCAC.a6.Type = Type(ac6_idx);
            GenSys.EXCAC.a6.GenNum = exc_con(ac6_idx,2);
            GenSys.EXCAC.a6.ExcIndex = ac6_idx;
            GenSys.EXCAC.a6.Ka = exc_con(ac6_idx,6);
            GenSys.EXCAC.a6.Ta = exc_con(ac6_idx,7);
            GenSys.EXCAC.a6.Tb = exc_con(ac6_idx,8);
            GenSys.EXCAC.a6.Tc = exc_con(ac6_idx,9);
            GenSys.EXCAC.a6.Tk = exc_con(ac6_idx,10);
            GenSys.EXCAC.a6.Tj = exc_con(ac6_idx,11);
            GenSys.EXCAC.a6.Kh = exc_con(ac6_idx,12);
            GenSys.EXCAC.a6.Th = exc_con(ac6_idx,13);
            GenSys.EXCAC.a6.VrMx = exc_con(ac6_idx,18);
            GenSys.EXCAC.a6.VrMn = exc_con(ac6_idx,19);
            GenSys.EXCAC.a6.VaMx = exc_con(ac6_idx,16);
            GenSys.EXCAC.a6.VaMn = exc_con(ac6_idx,17);
        end
        if~isempty(ac7_idx)
            GenNum = exc_con(ac7_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a7.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a7.PSSD = [];
            end
            GenSys.EXCAC.e7.Type = Type(ac7_idx);
            GenSys.EXCAC.e7.GenNum = exc_con(ac7_idx,2);
            GenSys.EXCAC.e7.ExcIndex = ac7_idx;
            GenSys.EXCAC.e7.Kc = exc_con(ac7_idx,20);
            GenSys.EXCAC.e7.Ke = exc_con(ac7_idx,21);
            GenSys.EXCAC.e7.Te = exc_con(ac7_idx,22);
            GenSys.EXCAC.e7.se1 = exc_con(ac7_idx,23);
            GenSys.EXCAC.e7.se2 = exc_con(ac7_idx,24);
            GenSys.EXCAC.e7.ve1 = exc_con(ac7_idx,25);
            GenSys.EXCAC.e7.ve2 = exc_con(ac7_idx,26);
            GenSys.EXCAC.e7.KL = exc_con(ac7_idx,27);
            GenSys.EXCAC.e7.Kp = exc_con(ac7_idx,28);
            GenSys.EXCAC.e7.Kd = exc_con(ac7_idx,29);
            GenSys.EXCAC.e7.Vfemx = exc_con(ac7_idx,30);                               
            GenSys.EXCAC.e7.VeMn = exc_con(ac7_idx,31);
            GenSys.EXCAC.a7.TransD.Type = Type(ac7_idx);
            GenSys.EXCAC.a7.TransD.GenNum = exc_con(ac7_idx,2);
            GenSys.EXCAC.a7.TransD.ExcIndex = ac7_idx; 
            GenSys.EXCAC.a7.TransD.Rc = Rc(ac7_idx);
            GenSys.EXCAC.a7.TransD.Xc = Xc(ac7_idx);
            GenSys.EXCAC.a7.TransD.Tr = Tr(ac7_idx);
            GenSys.EXCAC.a7.Type = Type(ac7_idx);
            GenSys.EXCAC.a7.GenNum = exc_con(ac7_idx,2);
            GenSys.EXCAC.a7.ExcIndex = ac7_idx;
            GenSys.EXCAC.a7.VrMx = exc_con(ac7_idx,18);
            GenSys.EXCAC.a7.VrMn = exc_con(ac7_idx,19);
            GenSys.EXCAC.a7.VaMx = exc_con(ac7_idx,16);
            GenSys.EXCAC.a7.VaMn = exc_con(ac7_idx,17);
            GenSys.EXCAC.a7.Kpr = exc_con(ac7_idx,6);
            GenSys.EXCAC.a7.Kir = exc_con(ac7_idx,7);
            GenSys.EXCAC.a7.Kpa = exc_con(ac7_idx,8);
            GenSys.EXCAC.a7.Kia = exc_con(ac7_idx,9);
            GenSys.EXCAC.a7.Kda = exc_con(ac7_idx,10);
            GenSys.EXCAC.a7.Tda = exc_con(ac7_idx,11);
            GenSys.EXCAC.a7.Kf1 = exc_con(ac7_idx,12);
            GenSys.EXCAC.a7.Kf2 = exc_con(ac7_idx,13);
            GenSys.EXCAC.a7.Kf3 = exc_con(ac7_idx,14);
            GenSys.EXCAC.a7.Tf = exc_con(ac7_idx,15);
        end
        if~isempty(ac8_idx)
            GenNum = exc_con(ac8_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCAC.a8.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCAC.a8.PSSD = [];
            end
            GenSys.EXCAC.e8.Type = Type(ac8_idx);
            GenSys.EXCAC.e8.GenNum = exc_con(ac8_idx,2);
            GenSys.EXCAC.e8.ExcIndex = ac8_idx;
            GenSys.EXCAC.e8.Kc = exc_con(ac8_idx,20);
            GenSys.EXCAC.e8.Ke = exc_con(ac8_idx,21);
            GenSys.EXCAC.e8.Te = exc_con(ac8_idx,22);
            GenSys.EXCAC.e8.se1 = exc_con(ac8_idx,23);
            GenSys.EXCAC.e8.se2 = exc_con(ac8_idx,24);
            GenSys.EXCAC.e8.ve1 = exc_con(ac8_idx,25);
            GenSys.EXCAC.e8.ve2 = exc_con(ac8_idx,26);
            GenSys.EXCAC.e8.Kd = exc_con(ac8_idx,29);
            GenSys.EXCAC.e8.Vfemx = exc_con(ac8_idx,30);                               
            GenSys.EXCAC.e8.VeMn = exc_con(ac8_idx,31);
            GenSys.EXCAC.a8.TransD.Type = Type(ac8_idx);
            GenSys.EXCAC.a8.TransD.GenNum = exc_con(ac8_idx,2);
            GenSys.EXCAC.a8.TransD.ExcIndex = ac8_idx; 
            GenSys.EXCAC.a8.TransD.Rc = Rc(ac8_idx);
            GenSys.EXCAC.a8.TransD.Xc = Xc(ac8_idx);
            GenSys.EXCAC.a8.TransD.Tr = Tr(ac8_idx);
            GenSys.EXCAC.a8.Type = Type(ac8_idx);
            GenSys.EXCAC.a8.GenNum = exc_con(ac8_idx,2);
            GenSys.EXCAC.a8.ExcIndex = ac8_idx;
            GenSys.EXCAC.a8.VrMx = exc_con(ac8_idx,18);
            GenSys.EXCAC.a8.VrMn = exc_con(ac8_idx,19);
            GenSys.EXCAC.a8.Ka = exc_con(ac8_idx,6);
            GenSys.EXCAC.a8.Ta = exc_con(ac8_idx,7);
            GenSys.EXCAC.a8.Kpa = exc_con(ac8_idx,8);
            GenSys.EXCAC.a8.Kia = exc_con(ac8_idx,9);
            GenSys.EXCAC.a8.Kda = exc_con(ac8_idx,10);
            GenSys.EXCAC.a8.Tda = exc_con(ac8_idx,11);
            
        end	
        GenSys.EXCST.e1 = [];GenSys.EXCST.a1 = [];
        GenSys.EXCST.e2 = [];GenSys.EXCST.a2 = [];
        GenSys.EXCST.e3 = [];GenSys.EXCST.a3 = [];
        GenSys.EXCST.e4 = [];GenSys.EXCST.a4 = [];
        GenSys.EXCST.e5 = [];GenSys.EXCST.a5 = [];
        if ~isempty(st1_idx)
            GenSys.EXCST.e1.Type = Type(st1_idx);
            GenNum = exc_con(st1_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCST.a1.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCST.a1.PSSD = [];
            end
            GenSys.EXCST.e1.GenNum = GenNum;
            GenSys.EXCST.e1.ExcIndex = st1_idx;
            GenSys.EXCST.e1.Kc = exc_con(st1_idx,20);
            GenSys.EXCST.a1.TransD.Type = Type(st1_idx);
            GenSys.EXCST.a1.TransD.GenNum = GenNum;
            GenSys.EXCST.a1.TransD.ExcIndex = st1_idx; 
            GenSys.EXCST.a1.TransD.Rc = Rc(st1_idx);
            GenSys.EXCST.a1.TransD.Xc = Xc(st1_idx);
            GenSys.EXCST.a1.TransD.Tr = Tr(st1_idx);
            GenSys.EXCST.a1.Type = Type(st1_idx);
            GenSys.EXCST.a1.Ka = exc_con(st1_idx,6);
            GenSys.EXCST.a1.Ta = exc_con(st1_idx,7);
            GenSys.EXCST.a1.Tb = exc_con(st1_idx,8);
            GenSys.EXCST.a1.Tc = exc_con(st1_idx,9);
            GenSys.EXCST.a1.Kf = exc_con(st1_idx,10);
            GenSys.EXCST.a1.Tf = exc_con(st1_idx,11);
            GenSys.EXCST.a1.Tb1 = exc_con(st1_idx,12);
            GenSys.EXCST.a1.Tc1 = exc_con(st1_idx,13);
            GenSys.EXCST.a1.ViMx = exc_con(st1_idx,14);
            GenSys.EXCST.a1.ViMn = exc_con(st1_idx,15);
            GenSys.EXCST.a1.VrMx = exc_con(st1_idx,18);
            GenSys.EXCST.a1.VrMn = exc_con(st1_idx,19);
        end
        if ~isempty(st2_idx)
            GenNum = exc_con(st2_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCST.a2.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCST.a2.PSSD = [];
            end
            GenSys.EXCST.e2.cs.Type = Type(st2_idx);
            GenSys.EXCST.e2.cs.GenNum = GenNum;
            GenSys.EXCST.e2.cs.Kc = exc_con(st2_idx,20);
            GenSys.EXCST.e2.cs.Kp = exc_con(st2_idx,23);
            GenSys.EXCST.e2.cs.Ki = exc_con(st2_idx,24);
            GenSys.EXCST.e2.cs.thetap = exc_con(st2_idx,25);
            GenSys.EXCST.e2.cs.Xl = exc_con(st2_idx,26);
            GenSys.EXCST.e2.cs.Vbmx = exc_con(st2_idx,28);
            GenSys.EXCST.e2.Type = Type(st2_idx);
            GenSys.EXCST.e2.GenNum = GenNum;
            GenSys.EXCST.e2.ExcIndex = st2_idx;
            GenSys.EXCST.e2.EfdMx = exc_con(st2_idx,30);
            GenSys.EXCST.e2.EfdMn = zeros(length(st2_idx),1);
            GenSys.EXCST.e2.Ke = exc_con(st2_idx,21);
            GenSys.EXCST.e2.Te = exc_con(st2_idx,22);
            GenSys.EXCST.a2.TransD.Type = Type(st2_idx);
            GenSys.EXCST.a2.TransD.GenNum = GenNum;
            GenSys.EXCST.a2.TransD.ExcIndex = st2_idx; 
            GenSys.EXCST.a2.TransD.Rc = Rc(st2_idx);
            GenSys.EXCST.a2.TransD.Xc = Xc(st2_idx);
            GenSys.EXCST.a2.TransD.Tr = Tr(st2_idx);
            GenSys.EXCST.a2.Type = Type(st2_idx);
            GenSys.EXCST.a2.Ka = exc_con(st2_idx,6);
            GenSys.EXCST.a2.Ta = exc_con(st2_idx,7);
            GenSys.EXCST.a2.Kf = exc_con(st2_idx,10);
            GenSys.EXCST.a2.Tf = exc_con(st2_idx,11);
            GenSys.EXCST.a2.VrMx = sparse(exc_con(st2_idx,18));
            GenSys.EXCST.a2.VrMn = sparse(exc_con(st2_idx,19));
        end
        if ~isempty(st3_idx)
            GenNum = exc_con(st3_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCST.a3.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCST.a3.PSSD = [];
            end
            GenSys.EXCST.e3.cs.Type = Type(st3_idx);
            GenSys.EXCST.e3.cs.GenNum = GenNum;
            GenSys.EXCST.e3.cs.Kc = exc_con(st3_idx,20);
            GenSys.EXCST.e3.cs.Kp = exc_con(st3_idx,23);
            GenSys.EXCST.e3.cs.Ki = exc_con(st3_idx,24);
            GenSys.EXCST.e3.cs.thetap = exc_con(st3_idx,25);
            GenSys.EXCST.e3.cs.Xl = exc_con(st3_idx,26);
            GenSys.EXCST.e3.cs.Vbmx = exc_con(st3_idx,28);
            GenSys.EXCST.e3.Type = Type(st3_idx);
            GenSys.EXCST.e3.GenNum = GenNum;
            GenSys.EXCST.e3.ExcIndex = st3_idx;
            GenSys.EXCST.e3.Kg = exc_con(st3_idx,27);
            GenSys.EXCST.e3.Vgmx = exc_con(st3_idx,29);
            GenSys.EXCST.e3.Ke = exc_con(st3_idx,21);
            GenSys.EXCST.e3.Te = exc_con(st3_idx,22);
            GenSys.EXCST.e3.Vemx = exc_con(st3_idx,30);
            GenSys.EXCST.e3.Vemn = exc_con(st3_idx,31); 
            GenSys.EXCST.a3.TransD.Type = Type(st3_idx);
            GenSys.EXCST.a3.TransD.GenNum = GenNum;
            GenSys.EXCST.a3.TransD.ExcIndex = st3_idx; 
            GenSys.EXCST.a3.TransD.Rc = Rc(st3_idx);
            GenSys.EXCST.a3.TransD.Xc = Xc(st3_idx);
            GenSys.EXCST.a3.TransD.Tr = Tr(st3_idx);
            GenSys.EXCST.a3.Type = Type(st3_idx);
            GenSys.EXCST.a3.Ka = exc_con(st3_idx,6);
            GenSys.EXCST.a3.Ta = exc_con(st3_idx,7);
            GenSys.EXCST.a3.Tb = exc_con(st3_idx,8);
            GenSys.EXCST.a3.Tc = exc_con(st3_idx,9);
            GenSys.EXCST.a3.ViMx = exc_con(st3_idx,14);
            GenSys.EXCST.a3.ViMn = exc_con(st3_idx,15);
            GenSys.EXCST.a3.VrMx = exc_con(st3_idx,18);
            GenSys.EXCST.a3.VrMn = exc_con(st3_idx,19);
        end
        if ~isempty(st4_idx)
            GenNum = exc_con(st4_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCST.a4.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCST.a4.PSSD = [];
            end
            GenSys.EXCST.e4.cs.Type = Type(st4_idx);
            GenSys.EXCST.e4.cs.GenNum = GenNum;
            GenSys.EXCST.e4.cs.Kc = exc_con(st4_idx,20);
            GenSys.EXCST.e4.cs.Kp = exc_con(st4_idx,23);
            GenSys.EXCST.e4.cs.Ki = exc_con(st4_idx,24);
            GenSys.EXCST.e4.cs.thetap = exc_con(st4_idx,25);
            GenSys.EXCST.e4.cs.Xl = exc_con(st4_idx,26);
            GenSys.EXCST.e4.cs.Vbmx = exc_con(st4_idx,28);
            GenSys.EXCST.e4.Type = Type(st4_idx);
            GenSys.EXCST.e4.GenNum = GenNum;
            GenSys.EXCST.e4.ExcIndex = st4_idx;
            GenSys.EXCST.e4.Kg = exc_con(st4_idx,27);
            GenSys.EXCST.e4.Kpe = exc_con(st4_idx,21);
            GenSys.EXCST.e4.Kie = exc_con(st4_idx,22);
            GenSys.EXCST.e4.Vemx = exc_con(st4_idx,30);
            GenSys.EXCST.e4.Vemn = exc_con(st4_idx,31); 
            GenSys.EXCST.a4.TransD.Type = Type(st4_idx);
            GenSys.EXCST.a4.TransD.GenNum = GenNum;
            GenSys.EXCST.a4.TransD.ExcIndex = st4_idx; 
            GenSys.EXCST.a4.TransD.Rc = Rc(st4_idx);
            GenSys.EXCST.a4.TransD.Xc = Xc(st4_idx);
            GenSys.EXCST.a4.TransD.Tr = Tr(st4_idx);
            GenSys.EXCST.a4.Type = Type(st4_idx);
            GenSys.EXCST.a4.Kpr = exc_con(st4_idx,6);
            GenSys.EXCST.a4.Ta = exc_con(st4_idx,8);
            GenSys.EXCST.a4.Kir=exc_con(st4_idx,7);
            GenSys.EXCST.a4.VrMx = exc_con(st4_idx,18);
            GenSys.EXCST.a4.VrMn = exc_con(st4_idx,19); 
        end
        if ~isempty(st5_idx)
            GenNum = exc_con(st5_idx,2);
            pss_idx = com_index(PGenNum,GenNum);
            if ~isempty(pss_idx)
                GenSys.EXCST.a5.PSSD = pdata(pss_con,pss_idx);
            else 
                GenSys.EXCST.a5.PSSD = [];
            end
            GenSys.EXCST.e5.Type = Type(st5_idx);
            GenSys.EXCST.e5.GenNum = GenNum(st5_idx);
            GenSys.EXCST.e5.ExcIndex = st5_idx;
            GenSys.EXCST.e5.Kc = exc_con(st5_idx,20);
            GenSys.EXCST.a5.TransD.Type = Type(st5_idx);
            GenSys.EXCST.a5.TransD.GenNum = exc_con(st5_idx,2);
            GenSys.EXCST.a5.TransD.ExcIndex = st5_idx; 
            GenSys.EXCST.a5.TransD.Rc = Rc(st5_idx);
            GenSys.EXCST.a5.TransD.Xc = Xc(st5_idx);
            GenSys.EXCST.a5.TransD.Tr = Tr(st5_idx);
            GenSys.EXCST.a5.Type = Type(st5_idx);
            GenSys.EXCST.a5.ExcIndex = st5_idx;
            GenSys.EXCST.a5.Ka = exc_con(st5_idx,6);
            GenSys.EXCST.a5.Ta = exc_con(st5_idx,7);
            GenSys.EXCST.a5.Tb = exc_con(st5_idx,8);
            GenSys.EXCST.a5.Tc = exc_con(st5_idx,9);
            GenSys.EXCST.a5.Tb1 = exc_con(st5_idx,12);
            GenSys.EXCST.a5.Tc1 = exc_con(st5_idx,13);
            GenSys.EXCST.a5.VrMx = exc_con(st5_idx,18);
            GenSys.EXCST.a5.VrMn = exc_con(st5_idx,19);
        end     
    end
    if ~isempty(hgt_con)
        n_hg = length(hgt_con(:,1));
        GenSys.HGovD.GenNum = hgt_con(:,1);
        GenSys.HGovD.Rp = hgt_con(:,2);
        GenSys.HGovD.Rt = hgt_con(:,3);
        GenSys.HGovD.Tr = hgt_con(:,4);
        GenSys.HGovD.Tp = hgt_con(:,5);
        GenSys.HGovD.Ks = hgt_con(:,6);
        GenSys.HGovD.Tg = hgt_con(:,7);
        GenSys.HGovD.MaxOr = hgt_con(:,8);
        GenSys.HGovD.MaxCr = hgt_con(:,9);
        GenSys.HGovD.db = hgt_con(:,10); 
        GenSys.HGovD.Rating = hgt_con(:,11);
        GenSys.HGovD.TW = hgt_con(:,12);
        GenSys.HGovD.gnl = hgt_con(:,13);
        GenSys.HGovD.gfl = hgt_con(:,14);
        GenSys.HGovD.H = hgt_con(:,15);
    else
        n_hg=0;
    end
    if ~isempty(tgt_con)
        if n_hg~=0
            ctb_idx = com_index(tgt_con(:,1),hgt_con(:,1));
            if ~isempty(ctb_idx)
                uiwait(msgbox('both thermal and hydro turbines specified for some generators','governor data error','modal'))
                tgt_con(ctb_idx,:)=[];
            end
        end
        GenSys.TGovD.GenNum = tgt_con(:,1);
        GenSys.TGovD.droop = tgt_con(:,2);
        GenSys.TGovD.Ttd = tgt_con(:,3);
        GenSys.TGovD.Tcs = tgt_con(:,4);
        GenSys.TGovD.Rcmax = tgt_con(:,5);
        GenSys.TGovD.Rcmin = tgt_con(:,6);
        GenSys.TGovD.ivdroop = tgt_con(:,7);
        GenSys.TGovD.ivref = tgt_con(:,8);
        GenSys.TGovD.Tis = tgt_con(:,9);
        GenSys.TGovD.Rimax = tgt_con(:,10);
        GenSys.TGovD.Rimin = tgt_con(:,11);
        GenSys.TGovD.Rating = tgt_con(:,12);
        GenSys.TGovD.THP = tgt_con(:,13);
        GenSys.TGovD.TRH = tgt_con(:,14);
        GenSys.TGovD.TLP1 = tgt_con(:,15);
        GenSys.TGovD.TLP2 = tgt_con(:,16); 
        GenSys.TGovD.FHP = tgt_con(:,17);
        GenSys.TGovD.FLP1 = tgt_con(:,18);
        GenSys.TGovD.FLP2 = tgt_con(:,19);
    end
end	

if (isempty(dcl_con)||isempty(dcr_con)||isempty(dci_con)||isempty(dcrc_con)||isempty(dcic_con));
    HVDCCD = [];
else
	HVDCCD.Line.LineNum = dcl_con(:,1);
	HVDCCD.Line.Ldcr = dcld_con(:,1);
	HVDCCD.Line.Ldci = dcld_con(:,2);
	HVDCCD.Line.Cl   = dcld_con(:,3);
	HVDCCD.Rectifier.Type = dcrc_con(:,1);
	HVDCCD.Rectifier.Tri = dcrc_con(:,2);
	HVDCCD.Rectifier.Trv = dcrc_con(:,3);
	HVDCCD.Rectifier.Ki = dcrc_con(:,4);
	HVDCCD.Rectifier.Kp = dcrc_con(:,5);
	HVDCCD.Rectifier.Ko = dcrc_con(:,6);
	HVDCCD.Rectifier.alphamax = dcrc_con(:,7);
	HVDCCD.Rectifier.alphamin = dcrc_con(:,8);
	HVDCCD.Inverter.Type = dcic_con(:,1);
	HVDCCD.Inverter.Tri = dcic_con(:,2);
	HVDCCD.Inverter.Trv = dcic_con(:,3);
	HVDCCD.Inverter.Ki = dcic_con(:,4);
	HVDCCD.Inverter.Kp = dcic_con(:,5);
	HVDCCD.Inverter.Ko = dcic_con(:,6);
	HVDCCD.Inverter.gammamax = dcic_con(:,7);
	HVDCCD.Inverter.gammamin = dcic_con(:,8);
end	
	
if ~isempty(svc_con)
    SVCD.BusNum = svc_con(:,2);
    SVCD.BasMva = svc_con(:,3);
    SVCD.BMax = svc_con(:,4);
    SVCD.BMin = svc_con(:,5);
    SVCD.Kr = svc_con(:,6);
    SVCD.Tr = svc_con(:,7);
    SVCD.Tb = svc_con(:,8);
    SVCD.Tc = svc_con(:,9); 
end
if ~isempty(tcsc_con)
    TCSCD.FromBusNum = tcsc_con(:,1);
    TCSCD.ToBusNum = tcsc_con(:,2);
    TCSCD.BasMva = tcsc_con(:,3);
    TCSCD.BMax = tcsc_con(:,4);
    TCSCD.BMin = tcsc_con(:,5);
    TCSCD.Kr = tcsc_con(:,6);
    TCSCD.Tr = tcsc_con(:,7);
    TCSCD.Tb = tcsc_con(:,8);
    TCSCD.Tc = tcsc_con(:,9);
end
if ~isempty(upfc_con)
	UPFCD.Num = upfc_con(:,1);
	UPFCD.KrSep = upfc_con(:,15);
	UPFCD.KiSep = upfc_con(:,16);
	UPFCD.TrSep = upfc_con(:,17);
	UPFCD.KrSeq = upfc_con(:,18);
	UPFCD.KiSeq = upfc_con(:,19);
	UPFCD.TrSeq = upfc_con(:,20);
	UPFCD.KrShv = upfc_con(:,21);
	UPFCD.KiShv = upfc_con(:,22);
	UPFCD.TrShv = upfc_con(:,23);
	UPFCD.KrCapv = upfc_con(:,24);
	UPFCD.KiCapv = upfc_con(:,25);
	UPFCD.TrCapv = upfc_con(:,26);
end
    
return
function p = pdata(pss_con,p_idx)
n_pss = length(pss_con(p_idx,1));
[gnu,ign]=unique(pss_con(p_idx,2));
if length(ign)~=n_pss
    uiwait(msgbox('there are multiple pss at some generators','pss data error','modal'))
    p_idx=p_idx(ign);
end
Type = pss_con(p_idx,1);
GenNum = pss_con(p_idx,2);
p.Type = Type;
p.GenNum = GenNum;
p.Gain = pss_con(p_idx,3);
p.Tw = pss_con(p_idx,4);
p.Tn1 = pss_con(p_idx,5);
p.Td1 = pss_con(p_idx,6);
p.Tn2 = pss_con(p_idx,7);
p.Td2 = pss_con(p_idx,8);
p.MaxOut = pss_con(p_idx,9);
p.MinOut = pss_con(p_idx,10);
dpwf_idx = find(Type>=4);
p.dpwf_idx = dpwf_idx;
p.dpwf.Type = Type(dpwf_idx)-3;
p.dpwf.GenNum = GenNum(dpwf_idx);
p.dpwf.Twd = pss_con(p_idx(dpwf_idx),11);
p.dpwf.Tnf = pss_con(p_idx(dpwf_idx),12);
p.dpwf.Tdf = pss_con(p_idx(dpwf_idx),13);
p.dpwf.nnf = pss_con(p_idx(dpwf_idx),14);
p.dpwf.ndf = pss_con(p_idx(dpwf_idx),15);
p.dpwf.Gp = pss_con(p_idx(dpwf_idx),16);
return

