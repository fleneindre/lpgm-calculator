# lpgm-calculator


このクラスは、気象庁の長周期地震動階級をリアルタイムの生加速度データから計算するように設計されています。

また、水平成分の絶対速度応答スペクトルを取得することもできます。これは、0.2秒刻みで1.6秒から7.8秒の範囲です。

使用した方法は、*Nigam, Navin C. and Jennings, Paul C. (1969) Calculation of 
response spectra from strong-motion earthquake records.*からのものです。

読み取り値は、加速度計が建物の1階、または地面に接続された最も水平な面に配置されている場合にのみ有効です。
サンプルレートは、時間の経過とともに一定である必要があります。
加速度入力はフィルタリングされていない必要があります（重力補正なし）


This class is designed to compute the Japan Meteorological Agency Long-Period Ground Motion class, from real-time raw acceleration data.
It is also possible to retrieve the Absolute Velocity Response Spectrum of horizontal components, spanning between 1.6s and 7.8s with 0.2s increments.

The method used is from *Nigam, Navin C. and Jennings, Paul C. (1969) Calculation of response spectra from strong-motion earthquake records.*

The readings are only valid if the accelerometer is placed on the first floor of any building or on a solid surface tied to the ground as level as 
possible. The sample rate has to be constant over time. The acceleration input has to be unfiltered (ie. without gravity compensation)

## Dependencies

The dependencies need to be installed by running the following command:
次のコマンドを実行して、依存関係をインストールする必要があります。 

```
pip install numpy scipy dvg-ringbuffer
```