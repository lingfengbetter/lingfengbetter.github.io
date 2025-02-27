+++
date = '2025-02-27T14:04:06+08:00'
draft = false

katex=true

title = '超声成像-----Null Subtraction Imaging'

+++

### 超声成像-----Null Subtraction Imaging

#### 摘要

Null Subtraction Imaging 是一种非线性图像处理技术（一种变迹方案），通过在相同 RF 数据上使用不同的接收窗，可以在保持低旁瓣水平的同时提高横向分辨率。将三种不同加窗函数生成的图像进行组合，来形成一个旁瓣水平低且横向分辨率有明显改善的图像。

==优点==：横向分辨率相较于矩形窗有高达30倍的改进，并且旁瓣水平也有较大程度的降低，计算效率高，变迹值与RF数据无关，可以预先计算。

==缺点==：CNR下降，需要调整直流偏移常数(通常为0.01～0.1)。

==适合的应用==：主要应用于检测和跟踪小的散射目标，检测无回声和强回声区域的能力有限。

#### 原理

线性变迹（调幅）与非线性调幅：传统的线性调幅虽能抑制旁瓣，但会牺牲主瓣宽度并降低横向分辨率。非线性调幅方法（如双/三调幅）通过优化波束形成算法，在保持主瓣宽度不变的情况下，将旁瓣水平降低9-10dB。

超声换能器或阵列的远场波束模式与其孔径函数和调幅函数乘积的傅里叶变换相关。例如，矩形函数的波束模式为sinc函数。基于这些傅里叶变换对，具有阵列零均值加权的调幅函数会在sinc函数的主瓣方向处产生零点。未调幅的矩形窗口（无调幅）的波束宽度受衍射限制保持恒定。然而，若比较矩形调幅与“倒置”（加负号反转过来）零点波束的宽度，对于零点波束而言主波束的衰减速度更快。这种零点波束的特性可以在保持低旁瓣的同时显著提高横向分辨率。

![NSI_1.jpg](/images/NSI_1.jpg)

![NSI_2.png](/images/NSI_2.png)

如果将zero-mean频谱反转过来就能获得较窄的主瓣。

变迹数值计算如下：
$$
A_{{\mathrm {R1},i}} = \begin{cases} 
1, & 1 \leq i < \dfrac {N}{2} \\\
-1, & \dfrac {N}{2} \leq i < N 
\end{cases} \tag{2}
$$
$$
A_{{\mathrm {R2},i}} = A_{{\mathrm {R1},i}} + c \tag{3}
$$

总的计算如下：
$$
E_{\mathrm {NSI}} = \frac {E_{\mathrm {DC1}} + E_{\mathrm {DC2}}}{2} - E_{\mathrm {ZM}}\tag{4}
$$
归一化：
$$
E_{\mathrm {NSI,dB}} = 20 \log \_{10}(E_{\mathrm {NSI}}) - \max \left[ 20 \log \_{10}(E_{\mathrm {NSI}}) \right]
\tag{5}
$$


```matlab
function ApodProfile=ApodGen(TransPara,BeamformPara,SizeRF)
    Depth=SizeRF(1);
    EleCount=SizeRF(2);
    ApodProfile = zeros((Depth-1)*BeamformPara.AxialInterp+1,EleCount,EleCount,4);
    [X,Z] = meshgrid((0:EleCount-1),(0:((Depth-1)*BeamformPara.AxialInterp)));
    dZ = BeamformPara.SoS/BeamformPara.SamplingFreq/BeamformPara.AxialInterp/2;
    Zd = Z*dZ+BeamformPara.InitDepth/TransPara.CenterFrequency*BeamformPara.SoS;
    for k = 1:EleCount
        aprlimit=min(k-1,EleCount-k)*2;
        aprsize=min(max(round(round(Zd(:,1)./BeamformPara.FNum./TransPara.Pitch)/2)*2,2),aprlimit);
        for q=1:(Depth-1)*BeamformPara.AxialInterp+1
            dasapod=ones(1,aprsize(q));
            zmlapod=[-ones(1,aprsize(q)/2),ones(1,aprsize(q)/2)];
            dclapod=zmlapod+BeamformPara.DCOffset;
            dcrapod=flip(dclapod,2);
            if(k<EleCount/2+1)
                ApodProfile(q,(k-aprsize(q)/2):(k+aprsize(q)/2-1),k,1)=dclapod;
                ApodProfile(q,(k-aprsize(q)/2):(k+aprsize(q)/2-1),k,2)=dcrapod;
                ApodProfile(q,(k-aprsize(q)/2):(k+aprsize(q)/2-1),k,3)=zmlapod;
                ApodProfile(q,(k-aprsize(q)/2):(k+aprsize(q)/2-1),k,4)=dasapod;
            else
                ApodProfile(q,(k-aprsize(q)/2+1):(k+aprsize(q)/2),k,1)=dclapod;
                ApodProfile(q,(k-aprsize(q)/2+1):(k+aprsize(q)/2),k,2)=dcrapod;
                ApodProfile(q,(k-aprsize(q)/2+1):(k+aprsize(q)/2),k,3)=zmlapod;
                ApodProfile(q,(k-aprsize(q)/2+1):(k+aprsize(q)/2),k,4)=dasapod;
            end
        end
    end
end
```

成像代码：

```matlab
NSILog=20*log10(mean(squeeze(abs(abs(NSI(:,:,1,:))+abs(NSI(:,:,2,:))-2*abs(NSI(:,:,3,:)))).^2,3));
NSINorm=NSILog-max(NSILog,[],'all');
```

#### 注意：在博客中写公式时“_”和“\\\”前需要加“\”，也即"\\\\\\"。

参考：

1. https://kiwamizamurai.github.io/posts/2022-03-06/
2. https://ieeexplore.ieee.org/document/8493541
