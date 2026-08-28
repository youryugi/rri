# RRI on iRIC

# Japanese

## RRI on iRICについて

RRI on iRICは、土木研究所ICHARMが開発したRRIモデルをiRICで使えるよう、RIverlinkの旭氏が中心となって、RRIモデルを改変したプログラムです。その後、土木研究所ICHARMが開発している降雨-土砂流出モデル等を統合したものです。
iRICのポストプロセッサを使って計算結果の表示を行うなど、便利に利用できるようになりました。

## 動作環境

iRICに同封されているものは，Intel Open APIを使ってコンパイルをしています．

同環境で使用したい場合は例えば下記のページなどで環境構築をお願いします．

https://qiita.com/Kazutake/items/a069f86d21ca43b6c153

## 改良版をiRICで動かす

RRIソルバーは make_RRI_ifx_static.bat でコンパイルできます．もちろん，そのほかの方法でコンパイルいただいても大丈夫だと思います．

コンパイルするとrri.exeが生成されるので，iRICのインストールフォルダにある，solvers/RRIにコピーすると更新ソルバーが使えます．

条件設定ファイルdefinition.xmlを変更した際も同様です．

## 著作権について

本プログラムの原著作物の著作権は、土木研究所（ICHARM）にあり、本プログラムは、RRIの二次著作物にあたります。本プログラムの利用、改変等に関する規約はRRIモデルに関するもの（https://www.pwri.go.jp/icharm/research/rri/rri_contract_j.html）と同一です。


# English

## About RRI on iRIC

RRI on iRIC is a modified version of the RRI model, originally developed by ICHARM (International Centre for Water Hazard and Risk Management) at PWRI (Public Works Research Institute). This adaptation was primarily led by Mr. Asahi of RIverlink to make the RRI model compatible with iRIC. Subsequently, it has been integrated with rainfall-sediment runoff models and other functionalities developed by ICHARM. It now offers enhanced usability, including the display of calculation results using iRIC's post-processor.

## Operating Environment

The version packaged with iRIC is compiled using the Intel Open API.
If you wish to use it in the same environment, please follow the instructions on a resource like the following page for environment setup:
https://qiita.com/Kazutake/items/a069f86d21ca43b6c153

## Running the Improved Version with iRIC

You can compile the RRI solver using make_RRI_ifx_static.bat. Of course, you may also compile it using other methods.
Compiling generates rri.exe. Copying this to the solvers/RRI directory within your iRIC installation folder will update the solver.
The same procedure applies when you modify the definition file definition.xml.

## Distribution Notes for the 2026-08-24 Build

This build keeps the visualization, labeling, and stability improvements from the development version, while disabling the site-specific inundation debug outputs used during model investigation. Ordinary runs do not create debug_inundation.csv, debug_sed_budget.csv, debug_inundation_spread.csv, debug_riv_shrink.csv, debug_progress.csv, or debug_heartbeat.txt.

The channel blockage extensions are included as optional experimental settings, but their defaults are disabled to preserve compatibility with general RRI use. These settings include bedload pass-through under low remaining channel depth, conveyance reduction by bed aggradation, river discharge limitation by remaining channel capacity, adjacent non-river-cell overtopping distribution, and pre-slope-routing river-slope exchange.

For the Egashira-type bedload formula, this build evaluates the Shields parameter for the bedload-driving term using the mean grain size dsm, rather than each grain-size fraction dsi. This change is intended to suppress unrealistically large fine-fraction bedload transport caused by the fraction-wise Shields scaling. Users comparing sediment transport results with older builds should note that bedload rates and bed changes can differ even when the optional blockage settings are disabled.

## Copyright

The copyright of the original work of this program belongs to the Public Works Research Institute (ICHARM), and this program is a secondary derivative work of RRI. The terms concerning the use, modification, etc., of this program are the same as those for the RRI model: 
https://www.pwri.go.jp/icharm/research/rri/rri_contract_e.html
